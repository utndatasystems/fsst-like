#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include "BenchmarkDriver.hpp"
// -------------------------------------------------------------------------------------
class StdFindEngine : public Engine {
public:
   StdFindEngine(std::string_view pattern)
       : pattern(pattern)
       , decode_buffer(128) {}

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (block.GetRow(row_idx).find(pattern) != std::string_view::npos) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);

         // Decode the row.
         uint32_t ideal_buffer_size = block.decoder.GetIdealBufferSize(compressed_text.size());
         if (ideal_buffer_size > decode_buffer.size()) {
            decode_buffer.resize(ideal_buffer_size);
         }
         uint32_t decoded_size = block.decoder.Decode(compressed_text, decode_buffer);

         // Match.
         std::string_view text(decode_buffer.data(), decoded_size);
         if (text.find(pattern) != std::string_view::npos) {
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
class StdFindEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (std::count(pattern.begin(), pattern.end(), '%') == 2 &&
          pattern.find('_') == std::string::npos &&
          pattern.starts_with('%') &&
          pattern.ends_with('%')) {
         return std::make_unique<StdFindEngine>(pattern.substr(1, pattern.size() - 2));
      }
      return nullptr;
   }

   std::string GetName() final { return "std_find"; }
};
// -------------------------------------------------------------------------------------
