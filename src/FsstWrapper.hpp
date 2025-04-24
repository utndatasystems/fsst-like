#pragma once
// -------------------------------------------------------------------------------------
#include <iostream>
#include <span>
#include "Utility.hpp"
#include "fsst.h"
// -------------------------------------------------------------------------------------
#define FSST_CORRUPT 32774747032022883 /* 7-byte number in little endian containing "corrupt" */
// -------------------------------------------------------------------------------------
class FsstEncoder : public NonCopyable {
public:
   ~FsstEncoder();

   // Interface to encode multiple strings at once.
   void InitializeEncoderOnSample(uint32_t row_count, const char** texts, uint64_t* lengths);
   uint32_t EncodeData(uint32_t row_count,
                       const char** texts, uint64_t* lengths,          // Input tuples
                       std::span<char> output,                         // Output buffer
                       char** encoded_texts, uint64_t* encoded_lengths // Output tuples
   ) const;

   // Write encode into flat output.
   static uint32_t GetRequiredDecoderSize();
   uint32_t SerializeDecoder(std::span<char> output) const;

   uint64_t GetEncodedSize(uint32_t encoded_row_count, char** encoded_texts, uint64_t* encoded_lengths) const;

private:
   fsst_encoder_t* encoder = nullptr; // Opaque pointer to fsst encoder.
};
// -------------------------------------------------------------------------------------
class FsstDecoder : public NonCopyable {
public:
   FsstDecoder() = default;
   FsstDecoder(FsstDecoder&& other);
   ~FsstDecoder();

   uint32_t DeserializeDecoder(std::span<const char> input);

   // Returns number of bytes written to output (failed if `>output.size()`).
   uint32_t Decode(std::span<const char> input, std::span<char> output) const;

   // We can also encode using the decoder.
   // After the encoder has been serialized, we only have the decoder object.
   // However, for equality restrictions, we need to encode the constant
   // to be able to compare with encoded tuple data. Hence, this method.
   std::pair<bool, uint32_t> Encode(std::string_view text, std::span<char> output) const;

   std::string SymbolToStr(unsigned code_index) const;
   void PrintSymbolTable(std::ostream& os) const;

   uint32_t GetSymbolTableSize() const { return symbol_table_size; }
   fsst_decoder_t& GetSymbolTable() const { return *decoder; }

   // How much memory is required so that FSST can fast decompress into it.
   uint32_t GetIdealBufferSize(uint32_t compressed_size) const { return compressed_size * 8 + 32; }

   uint8_t FindLongestSymbol(std::string_view text, bool allow_prefix = false) const;

   // Iterate a FSST-encoded string.
   // TODO: Make sure the latest bugs have been fixed from cwida/fsst.
   template<typename ConsumeCode, typename ConsumeChar>
   inline bool Iterate(
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

      // TODO: Remove the stuff with `posOut` (like everything).
      // TODO Also, _maybe_, the part with `size`. Not sure about the `if`-stmt for the last 3 bytes.
      while (posOut+32 <= size && posIn+4 <= lenIn) {
         unsigned int nextBlock, escapeMask;
         memcpy(&nextBlock, strIn+posIn, sizeof(unsigned int));
         escapeMask = (nextBlock&0x80808080u)&((((~nextBlock)&0x7F7F7F7Fu)+0x7F7F7F7Fu)^0x80808080u);
         if (escapeMask == 0) {
            code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
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
               // Consume the regular byte.
               if (!consume_char(strIn[posIn + 1])) return true;
               
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
         assert(0); // strOut[size-1] = 0;
      }

      return false; /* full size of decompressed string (could be >size, then the actually decompressed part) */
   }
private:
   uint32_t symbol_table_size;
   fsst_decoder_t* decoder = nullptr;

   static uint32_t CountMatchingBytes(std::string_view text, uint64_t symbol, uint32_t len);
};
// -------------------------------------------------------------------------------------
