#pragma once
// -------------------------------------------------------------------------------------
#include "BenchmarkDriver.hpp"
#include "TokenizerCodec.hpp"
// -------------------------------------------------------------------------------------
class EqualsEngine : public Engine {
public:
   explicit EqualsEngine(std::string_view pattern)
       : pattern(pattern)
       , decode_buffer(128)
       , encoded_pattern_buffer(128)
   {
      std::span<const char> pattern_span(pattern.data(), pattern.size());
      tokenizer_codec::Compress(pattern_span, tokenizer_encoded_pattern);
   }

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) final
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (block.GetRow(row_idx) == pattern) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const CompressedBlock& block, std::vector<uint32_t>& result) final
   {
      uint32_t match_count = 0;
      bool has_encoded_pattern = false;
      uint32_t encoded_pattern_size = 0;

      if (!block.decoder.IsTokenizerMode()) {
         const uint32_t required = static_cast<uint32_t>(pattern.size() * 2 + 8);
         if (required > encoded_pattern_buffer.size()) {
            encoded_pattern_buffer.resize(required);
         }
         auto [ok, written] = block.decoder.Encode(pattern, encoded_pattern_buffer);
         has_encoded_pattern = ok;
         encoded_pattern_size = written;
         assert(has_encoded_pattern && "FSST equality fast path expects full predicate encoding");
      }

      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         std::string_view compressed_text = block.GetRow(row_idx);
         if (CompressedEquals(block.decoder, compressed_text, has_encoded_pattern, encoded_pattern_size)) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

private:
   bool CompressedEquals(const FsstDecoder& decoder, std::string_view compressed_text, bool has_encoded_pattern, uint32_t encoded_pattern_size)
   {
      if (decoder.IsTokenizerMode()) {
         if (compressed_text.size() != tokenizer_encoded_pattern.size()) {
            return false;
         }
         return std::memcmp(compressed_text.data(), tokenizer_encoded_pattern.data(), tokenizer_encoded_pattern.size()) == 0;
      }

      if (has_encoded_pattern && compressed_text.size() == encoded_pattern_size) {
         return std::memcmp(compressed_text.data(), encoded_pattern_buffer.data(), encoded_pattern_size) == 0;
      }

      // Fallback: if the pattern cannot be encoded exactly with this symbol table,
      // decode and compare in text space.
      uint32_t ideal_buffer_size = decoder.GetIdealBufferSize(compressed_text.size());
      if (ideal_buffer_size > decode_buffer.size()) {
         decode_buffer.resize(ideal_buffer_size);
      }
      uint32_t decoded_size = decoder.Decode(compressed_text, decode_buffer);
      return std::string_view(decode_buffer.data(), decoded_size) == pattern;
   }

   std::string_view pattern;
   std::vector<char> decode_buffer;
   std::vector<char> encoded_pattern_buffer;
   std::vector<char> tokenizer_encoded_pattern;
};
// -------------------------------------------------------------------------------------
class EqualsEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (pattern.find('%') != std::string::npos || pattern.find('_') != std::string::npos) {
         return nullptr;
      }
      return std::make_unique<EqualsEngine>(pattern);
   }

   std::string GetName() final { return "equals"; }
};
// -------------------------------------------------------------------------------------
