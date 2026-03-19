#pragma once
// -------------------------------------------------------------------------------------
#include <cstring>
#include "BenchmarkDriver.hpp"
#include "codecs/TokenizerCodec.hpp"
// -------------------------------------------------------------------------------------
class EqualsEngine : public Engine {
public:
   explicit EqualsEngine(std::string_view pattern)
       : pattern(pattern)
       , encoded_pattern_buffer(128)
   {
      std::span<const char> pattern_span(pattern.data(), pattern.size());
      tokenizer_codec::Compress(pattern_span, tokenizer_encoded_pattern);
   }

   ~EqualsEngine() override
   {
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
      bool has_fsst_encoded_pattern = false;
      uint32_t fsst_encoded_pattern_size = 0;

      if (block.IsFsstCodec()) {
         const uint32_t required = static_cast<uint32_t>(pattern.size() * 2 + 8);
         if (required > encoded_pattern_buffer.size()) {
            encoded_pattern_buffer.resize(required);
         }
      auto [ok, written] = block.fsst_decoder.Encode(pattern, encoded_pattern_buffer);
      has_fsst_encoded_pattern = ok;
      fsst_encoded_pattern_size = written;
   }

      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         std::string_view compressed_text = block.GetRow(row_idx);
         if (CompressedEquals(block, compressed_text, has_fsst_encoded_pattern, fsst_encoded_pattern_size)) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

private:
   bool CompressedEquals(const CompressedBlock& block,
                         std::string_view compressed_text,
                         bool has_fsst_encoded_pattern,
                         uint32_t fsst_encoded_pattern_size)
   {
      if (block.IsTokenizerCodec()) {
         std::string_view encoded_pattern(tokenizer_encoded_pattern.data(), tokenizer_encoded_pattern.size());
         return compressed_text == encoded_pattern;
      }

      if (block.IsFsstCodec() == false) {
         const uint32_t ideal_buffer_size = block.GetIdealBufferSize(compressed_text.size());
         if (ideal_buffer_size > decode_buffer.size()) {
            decode_buffer.resize(ideal_buffer_size);
         }
         const uint32_t decoded_size = block.Decode(compressed_text, decode_buffer);
         std::string_view decoded_text(decode_buffer.data(), decoded_size);
         return decoded_text == pattern;
      }

      assert(has_fsst_encoded_pattern && "FSST equals expects pattern to always encode into pre-sized buffer");
      std::string_view encoded_pattern(encoded_pattern_buffer.data(), fsst_encoded_pattern_size);
      return compressed_text == encoded_pattern;
   }

   std::string_view pattern;
   std::vector<char> encoded_pattern_buffer;
   std::vector<char> tokenizer_encoded_pattern;
   std::vector<char> decode_buffer;
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
