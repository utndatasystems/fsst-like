#include "FsstWrapper.hpp"
#include "fsst.h"
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
static_assert(sizeof(uint64_t) == sizeof(size_t)); // FSST uses size_t everywhere.
// -------------------------------------------------------------------------------------
FsstEncoder::~FsstEncoder()
{
   if (encoder) {
      fsst_destroy(encoder);
      encoder = nullptr;
   }
}
// -------------------------------------------------------------------------------------
void FsstEncoder::InitializeEncoderOnSample(uint32_t row_count, const char** texts, uint64_t* lengths)
{
   assert(!encoder);

   // Cast inputs to FSST types.
   const size_t* length_array = reinterpret_cast<const size_t*>(lengths);
   const unsigned char** cast_texts = reinterpret_cast<const unsigned char**>(texts);

   // Sample (created encoder).
   encoder = fsst_create(row_count,
                         length_array,
                         cast_texts,
                         0);
}
// -------------------------------------------------------------------------------------
uint32_t FsstEncoder::EncodeData(uint32_t row_count, const char** texts, uint64_t* lengths, span<char> output, char** encoded_texts, uint64_t* encoded_lengths) const
{
   assert(encoder);

   // Cast inputs to FSST types.
   static_assert(sizeof(uint64_t) == sizeof(size_t));
   const size_t* length_array = reinterpret_cast<const size_t*>(lengths);
   size_t* encoded_length_array = reinterpret_cast<size_t*>(encoded_lengths);
   const unsigned char** cast_texts = reinterpret_cast<const unsigned char**>(texts);
   unsigned char** cast_encoded_texts = reinterpret_cast<unsigned char**>(encoded_texts);
   unsigned char* cast_output = reinterpret_cast<unsigned char*>(output.data());

   // Encode.
   return fsst_compress(encoder,
                        row_count,
                        length_array,
                        cast_texts,
                        output.size(),
                        cast_output,
                        encoded_length_array,
                        cast_encoded_texts);
}
// -------------------------------------------------------------------------------------
uint64_t FsstEncoder::GetEncodedSize(uint32_t encoded_row_count, char** encoded_texts, uint64_t* encoded_lengths) const
{
   if (encoded_row_count == 0) {
      return 0;
   }
   return (encoded_texts[encoded_row_count - 1] - encoded_texts[0]) + encoded_lengths[encoded_row_count - 1];
}
// -------------------------------------------------------------------------------------
uint32_t FsstEncoder::GetRequiredDecoderSize()
{
   return FSST_MAXHEADER;
}
// -------------------------------------------------------------------------------------
uint32_t FsstEncoder::SerializeDecoder(span<char> output) const
{
   assert(encoder);
   assert(output.size() >= FSST_MAXHEADER);

   unsigned char* cast_output = reinterpret_cast<unsigned char*>(output.data());

   return fsst_export(encoder, cast_output);
}
// -------------------------------------------------------------------------------------
FsstDecoder::FsstDecoder(FsstDecoder&& other)
{
   if (decoder) {
      delete reinterpret_cast<fsst_decoder_t*>(decoder);
   }
   decoder = other.decoder;
   symbol_table_size = other.symbol_table_size;
   other.decoder = nullptr;
}
// -------------------------------------------------------------------------------------
FsstDecoder::~FsstDecoder()
{
   if (decoder) {
      delete reinterpret_cast<fsst_decoder_t*>(decoder);
      decoder = nullptr;
   }
}
// -------------------------------------------------------------------------------------
uint32_t FsstDecoder::DeserializeDecoder(span<const char> input)
{
   if (!decoder) {
      decoder = new fsst_decoder_t();
   }

   const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(input.data());

   size_t symbol_table_size = 0;
   uint32_t consumed = fsst_import(reinterpret_cast<fsst_decoder_t*>(decoder), cast_input, &symbol_table_size);
   this->symbol_table_size = symbol_table_size;
   assert(consumed > 0);

   return consumed;
}
// -------------------------------------------------------------------------------------
uint32_t FsstDecoder::Decode(span<const char> input, span<char> decoded) const
{
   assert(decoder);
   const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(input.data());
   unsigned char* cast_decoded = reinterpret_cast<unsigned char*>(decoded.data());

   return fsst_decompress(reinterpret_cast<const fsst_decoder_t*>(decoder),
                          input.size(),
                          cast_input,
                          decoded.size(),
                          cast_decoded);
}
// -------------------------------------------------------------------------------------
pair<bool, uint32_t> FsstDecoder::Encode(string_view text, span<char> output) const
{
   fsst_decoder_t* symbol_table = reinterpret_cast<fsst_decoder_t*>(decoder);

   uint32_t read_idx = 0;
   uint32_t write_idx = 0;

   while (read_idx < text.size() && write_idx + 1 < output.size()) {
      uint32_t longest_symbol = FindLongestSymbol(text.substr(read_idx));
      if (longest_symbol >= 255) {
         output[write_idx++] = static_cast<char>(255);
         output[write_idx++] = text[read_idx++];
      }
      else {
         output[write_idx++] = static_cast<char>(longest_symbol);
         read_idx += symbol_table->len[longest_symbol];
      }
   }

   return make_pair(read_idx == text.size(), write_idx);
}
// -------------------------------------------------------------------------------------
void FsstDecoder::PrintSymbolTable(ostream& os) const
{
   fsst_decoder_t* symbol_table = reinterpret_cast<fsst_decoder_t*>(decoder);
   for (uint32_t idx = 0; idx < symbol_table_size; idx++) {
      os << "idx: " << idx << ", len: " << static_cast<int>(symbol_table->len[idx]) << ", symbol: ";
      for (uint32_t jdx = 0; jdx < symbol_table->len[idx]; jdx++) {
         os << static_cast<char>(symbol_table->symbol[idx] >> (jdx * 8));
      }
      os << endl;
   }
}
// -------------------------------------------------------------------------------------
uint8_t FsstDecoder::FindLongestSymbol(string_view text, bool allow_prefix) const
{
   fsst_decoder_t* symbol_table = reinterpret_cast<fsst_decoder_t*>(decoder);

   uint32_t longest_match = 0;
   uint32_t longest_match_idx = 255;
   for (uint32_t idx = 0; idx < symbol_table_size; idx++) {
      uint32_t matching_bytes = CountMatchingBytes(text, symbol_table->symbol[idx], symbol_table->len[idx]);
      if (allow_prefix || matching_bytes == symbol_table->len[idx]) {
         if (matching_bytes > longest_match) {
            longest_match = matching_bytes;
            longest_match_idx = idx;
         }
      }
   }
   if (longest_match == 0) {
      return 255;
   }
   return longest_match_idx;
}
// -------------------------------------------------------------------------------------
uint32_t FsstDecoder::CountMatchingBytes(string_view text, uint64_t symbol, uint32_t len)
{
   const char* symbol_str = reinterpret_cast<const char*>(&symbol);
   uint32_t count = 0;
   for (uint32_t idx = 0; idx < len && idx < text.size(); idx++) {
      if (symbol_str[idx] == text[idx]) {
         count++;
      }
      else {
         return count;
      }
   }
   return count;
}
// -------------------------------------------------------------------------------------
template<typename ConsumeCode, typename ConsumeChar>
inline bool FsstDecoder::Iterate(
  // IN: The byte-length of compressed string.
  size_t lenIn,
  // IN: The compressed string.
  const unsigned char *strIn,
  // IN: The byte-length of output buffer.
  size_t size,
  // IN: The char consumer.
  ConsumeChar&& consume_char,
  // IN: The code consumer.
  ConsumeCode&& consume_code
) const {
   unsigned char*__restrict__ len = (unsigned char* __restrict__) decoder->len;
  //  unsigned char*__restrict__ strOut = (unsigned char* __restrict__) output;
   unsigned long long*__restrict__ symbol = (unsigned long long* __restrict__) decoder->symbol; 
   size_t code, posOut = 0, posIn = 0;
#ifndef FSST_MUST_ALIGN /* defining on platforms that require aligned memory access may help their performance */
#define FSST_UNALIGNED_STORE(dst,src) memcpy((unsigned long long*) (dst), &(src), sizeof(unsigned long long))
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)

  // TODO: Remove the stuff with `posOut`. Like everything.
  // And even with `size` since we don't have any constraint right now.
   while (posOut+32 <= size && posIn+4 <= lenIn) {
      unsigned int nextBlock, escapeMask;
      memcpy(&nextBlock, strIn+posIn, sizeof(unsigned int));
      escapeMask = (nextBlock&0x80808080u)&((((~nextBlock)&0x7F7F7F7Fu)+0x7F7F7F7Fu)^0x80808080u);
      if (escapeMask == 0) {
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; //FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; //FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
     } else { 
         unsigned long firstEscapePos=__builtin_ctzl((unsigned long long) escapeMask)>>3;
         switch(firstEscapePos) { /* Duff's device */
         case 3: code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code];
                 // fall through
         case 2: code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code];
                 // fall through
         case 1: code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code];
                 // fall through
         case 0: posIn+=2; if (!consume_char(strIn[posIn - 1])) return true; posOut++; // strOut[posOut++] =  strIn[posIn-1]; /* decompress an escaped byte */
         }
      }
   }
   if (posOut+32 <= size) { // handle the possibly 3 last bytes without a loop
      if (posIn+2 <= lenIn) {
        if (strIn[posIn] != FSST_ESC) {
            code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            if (strIn[posIn] != FSST_ESC) {
               code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            } else { 
               posIn += 2; if (!consume_char(strIn[posIn - 1])) return true; posOut++; // strOut[posOut++] = strIn[posIn-1]; 
            }
         } else {
            // Consume the regualr byte.
            if (!consume_char(strIn[posIn + 1]))
              return true;
            
            posIn += 2; posOut++; // posOut++;
         } 
      }
      if (posIn < lenIn) { // last code cannot be an escape
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
      }
   }
#else
   while (posOut+8 <= size && posIn < lenIn)
      if ((code = strIn[posIn++]) < FSST_ESC) { /* symbol compressed as code? */
         FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); /* unaligned memory write */
         posOut += len[code];
      } else { 
         strOut[posOut] = strIn[posIn]; /* decompress an escaped byte */
         posIn++; posOut++; 
      }
#endif
#endif
   while (posIn < lenIn)
      if ((code = strIn[posIn++]) < FSST_ESC) {
         size_t posWrite = posOut, endWrite = posOut + len[code];
         unsigned char* __restrict__ symbolPointer = ((unsigned char* __restrict__) &symbol[code]) - posWrite;
         if ((posOut = endWrite) > size) endWrite = size;
         for(; posWrite < endWrite; posWrite++) { /* only write if there is room */
            if (!consume_char(symbolPointer[posWrite])) return true; // strOut[posWrite] = symbolPointer[posWrite];
         }
      } else {
         if (posOut < size) {
            if (!consume_char(strIn[posIn])) return true; // strOut[posOut] = strIn[posIn]; /* idem */
         }
         posIn++; posOut++;
      } 
   if (posOut >= size && (decoder->zeroTerminated&1)) {
     assert(0); //  strOut[size-1] = 0;
   }
  //  std::cerr << std::endl;
   return false; /* full size of decompressed string (could be >size, then the actually decompressed part) */
}