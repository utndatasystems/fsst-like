#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <string>
#include <vector>
#include "BenchmarkDriver.hpp"
// -------------------------------------------------------------------------------------
class OnPairLikeEngine : public Engine {
public:
   enum class Mode {
      Equals,
      Contains,
      StartsWith,
      Generic
   };

   OnPairLikeEngine(std::string_view pattern, Mode mode)
       : pattern(pattern)
       , mode(mode)
   {
      if (mode == Mode::Generic) {
         parts = SplitPattern(this->pattern);
      }
   }

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) final
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (Matches(block.GetRow(row_idx))) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const CompressedBlock& block, std::vector<uint32_t>& result) final
   {
      if (!block.IsOnPairCodec()) {
         return 0;
      }

      uint32_t match_count = 0;
      auto view = block.onpair_column.view();
      auto on_match = [&](size_t row_idx) {
         result[match_count++] = static_cast<uint32_t>(row_idx);
      };

      if (mode == Mode::Equals) {
         view.equals(pattern, on_match);
         return match_count;
      }
      if (mode == Mode::Contains) {
         view.contains(pattern, on_match);
         return match_count;
      }
      if (mode == Mode::StartsWith) {
         view.starts_with(pattern, on_match);
         return match_count;
      }

      const size_t buffer_size = static_cast<size_t>(block.max_uncompressed_row_size) + onpair::DECOMPRESS_BUFFER_PADDING;
      if (decode_buffer.size() < buffer_size) {
         decode_buffer.resize(buffer_size);
      }
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         const size_t decoded_size = view.decompress(row_idx, decode_buffer.data());
         if (Matches(std::string_view(decode_buffer.data(), decoded_size))) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

private:
   bool Matches(std::string_view text) const
   {
      if (mode == Mode::Equals) {
         return text == pattern;
      }
      if (mode == Mode::Contains) {
         return text.find(pattern) != std::string_view::npos;
      }
      if (mode == Mode::StartsWith) {
         return text.starts_with(pattern);
      }

      size_t start_pos = 0;
      if (!pattern.starts_with('%')) {
         if (parts.empty() || !text.starts_with(parts.front())) {
            return false;
         }
         start_pos = parts.front().size();
      }

      const size_t first_part = pattern.starts_with('%') ? 0 : 1;
      for (size_t idx = first_part; idx < parts.size(); idx++) {
         auto pos = text.find(parts[idx], start_pos);
         if (pos == std::string_view::npos) {
            return false;
         }
         start_pos = pos + parts[idx].size();
      }

      if (!pattern.ends_with('%') && !parts.empty()) {
         return text.ends_with(parts.back());
      }
      return true;
   }

   std::string pattern;
   Mode mode;
   std::vector<std::string_view> parts;
   std::vector<char> decode_buffer;
};
// -------------------------------------------------------------------------------------
class OnPairLikeEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (pattern.find('_') != std::string_view::npos) {
         return nullptr;
      }

      const size_t percent_count = std::count(pattern.begin(), pattern.end(), '%');
      if (percent_count == 0) {
         return std::make_unique<OnPairLikeEngine>(pattern, OnPairLikeEngine::Mode::Equals);
      }
      if (percent_count == 1 && pattern.ends_with('%')) {
         return std::make_unique<OnPairLikeEngine>(pattern.substr(0, pattern.size() - 1),
                                                   OnPairLikeEngine::Mode::StartsWith);
      }
      if (percent_count == 2 && pattern.starts_with('%') && pattern.ends_with('%')) {
         auto inner = pattern.substr(1, pattern.size() - 2);
         if (inner.find('%') == std::string_view::npos) {
            return std::make_unique<OnPairLikeEngine>(inner, OnPairLikeEngine::Mode::Contains);
         }
      }

      return std::make_unique<OnPairLikeEngine>(pattern, OnPairLikeEngine::Mode::Generic);
   }

   std::string GetName() final { return "onpair_like"; }
};
// -------------------------------------------------------------------------------------
