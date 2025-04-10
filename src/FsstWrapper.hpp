#pragma once
// -------------------------------------------------------------------------------------
#include <iostream>
#include <span>
#include "Utility.hpp"
#include "fsst.h"
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

   void PrintSymbolTable(std::ostream& os) const;

   uint32_t GetSymbolTableSize() const { return symbol_table_size; }
   fsst_decoder_t& GetSymbolTable() { return *decoder; }

   // How much memory is required so that FSST can fast decompress into it.
   uint32_t GetIdealBufferSize(uint32_t compressed_size) const { return compressed_size * 8 + 32; }

   uint8_t FindLongestSymbol(std::string_view text, bool allow_prefix = false) const;

private:
   uint32_t symbol_table_size;
   fsst_decoder_t* decoder = nullptr;

   static uint32_t CountMatchingBytes(std::string_view text, uint64_t symbol, uint32_t len);
};
// -------------------------------------------------------------------------------------
