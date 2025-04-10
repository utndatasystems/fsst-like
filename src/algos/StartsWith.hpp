#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include "BenchmarkDriver.hpp"
// -------------------------------------------------------------------------------------
class StartsWithEngine : public Engine {
public:
   StartsWithEngine(std::string_view pattern)
       : pattern(pattern)
       , decode_buffer(128)
   {
      assert(pattern.size() >= 9); // short ones not implemented .. there we might need to consider several symbols
   }

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (block.GetRow(row_idx).starts_with(pattern)) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result)
   {
      uint8_t symbol = block.decoder.FindLongestSymbol(pattern, false);

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);

         if (compressed_text.size() > 0 && symbol != compressed_text[0]) {
            continue;
         }

         // Decode the row.
         uint32_t ideal_buffer_size = block.decoder.GetIdealBufferSize(compressed_text.size());
         if (ideal_buffer_size > decode_buffer.size()) {
            decode_buffer.resize(ideal_buffer_size);
         }
         uint32_t decoded_size = block.decoder.Decode(compressed_text, decode_buffer);

         // Match.
         std::string_view text(decode_buffer.data(), decoded_size);
         if (text.starts_with(pattern)) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

private:
   std::string_view pattern;
   std::vector<char> decode_buffer;
};
// -------------------------------------------------------------------------------------
class StartsWithEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (std::count(pattern.begin(), pattern.end(), '%') == 1 &&
          pattern.find('_') == std::string::npos &&
          pattern.ends_with('%')) {
         return std::make_unique<StartsWithEngine>(pattern.substr(0, pattern.size() - 1));
      }
      return nullptr;
   }

   std::string GetName() final { return "starts_with"; }
};
// -------------------------------------------------------------------------------------
