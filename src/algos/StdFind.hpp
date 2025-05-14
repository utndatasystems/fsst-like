#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <ranges>
#include "BenchmarkDriver.hpp"
// -------------------------------------------------------------------------------------
template <class T>
class StdEngine : public Engine {
public:
   StdEngine(std::string_view pattern)
       : pattern(pattern)
       , decode_buffer(128) {}

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) final
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (static_cast<T*>(this)->Matches(block.GetRow(row_idx))) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result) final
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
         if (static_cast<T*>(this)->Matches(text)) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

protected:
   std::string_view pattern;
   std::vector<char> decode_buffer;
};
// -------------------------------------------------------------------------------------
class GeneralStdFindEngine : public StdEngine<GeneralStdFindEngine> {
public:
    explicit GeneralStdFindEngine(std::string_view pattern)
        : StdEngine(pattern) 
    {
        patterns_ = SplitPattern(pattern);
    }

    bool Matches(std::string_view text) const noexcept {
        std::size_t start_pos = 0;
        for (const auto& pat : patterns_) {
            auto pos = text.find(pat, start_pos);
            if (pos == std::string_view::npos)
                return false;

            // Advance after match.
            start_pos = pos + pat.size();
        }
        return true;
    }

private:
    std::vector<std::string_view> patterns_;
};
// -------------------------------------------------------------------------------------
class StdFindEngine : public StdEngine<StdFindEngine> {
public:
   StdFindEngine(std::string_view pattern)
       : StdEngine(pattern) {}

   bool Matches(std::string_view text) noexcept { return text.find(pattern) != std::string_view::npos; }
};
// -------------------------------------------------------------------------------------
class StdStartsWithEngine : public StdEngine<StdStartsWithEngine> {
public:
   StdStartsWithEngine(std::string_view pattern)
       : StdEngine(pattern) {}

   bool Matches(std::string_view text) noexcept { return text.starts_with(pattern); }
};
// -------------------------------------------------------------------------------------
class StdEndsWithEngine : public StdEngine<StdEndsWithEngine> {
public:
   StdEndsWithEngine(std::string_view pattern)
       : StdEngine(pattern) {}

   bool Matches(std::string_view text) noexcept { return text.ends_with(pattern); }
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
      if (std::count(pattern.begin(), pattern.end(), '%') == 1 &&
          pattern.find('_') == std::string::npos &&
          pattern.ends_with('%')) {
         return std::make_unique<StdStartsWithEngine>(pattern.substr(0, pattern.size() - 1));
      }
      if (std::count(pattern.begin(), pattern.end(), '%') == 1 &&
          pattern.find('_') == std::string::npos &&
          pattern.starts_with('%')) {
         return std::make_unique<StdEndsWithEngine>(pattern.substr(1, pattern.size() - 1));
      }

      // Generic case for: %p1%p2%.
      // TODO: We also need the more general case: [%]p1%p2[%].
      if (pattern.starts_with('%') &&
          pattern.ends_with('%') &&
          pattern.find('_') == std::string::npos) {
         return std::make_unique<GeneralStdFindEngine>(pattern.substr(1, pattern.size() - 2));
      }

      return nullptr;
   }

   std::string GetName() final { return "stl"; }
};
// -------------------------------------------------------------------------------------
