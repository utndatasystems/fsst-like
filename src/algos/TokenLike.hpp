#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>
#include "BenchmarkDriver.hpp"
#include "codecs/TokenizerDictionary.hpp"
// -------------------------------------------------------------------------------------
class PercentLikeMatcher {
public:
   explicit PercentLikeMatcher(std::string_view pattern)
       : pattern(pattern)
       , active(pattern.size() + 1)
       , next(pattern.size() + 1)
   {
   }

   void Reset()
   {
      std::fill(active.begin(), active.end(), false);
      active[0] = true;
      ApplyEpsilonClosure(active);
   }

   void Consume(std::string_view text)
   {
      for (char c : text) {
         Consume(c);
      }
   }

   bool Accepted() const
   {
      return active[pattern.size()];
   }

private:
   void Consume(char c)
   {
      std::fill(next.begin(), next.end(), false);
      for (size_t state = 0; state < pattern.size(); state++) {
         if (!active[state]) {
            continue;
         }
         if (pattern[state] == '%') {
            next[state] = true;
         } else if (pattern[state] == c) {
            next[state + 1] = true;
         }
      }

      ApplyEpsilonClosure(next);
      active.swap(next);
   }

   void ApplyEpsilonClosure(std::vector<bool>& states) const
   {
      for (size_t state = 0; state < pattern.size(); state++) {
         if (states[state] && pattern[state] == '%') {
            states[state + 1] = true;
         }
      }
   }

   std::string pattern;
   std::vector<bool> active;
   std::vector<bool> next;
};
// -------------------------------------------------------------------------------------
class TokenContainsMatcher {
public:
   explicit TokenContainsMatcher(std::string_view needle)
       : needle(needle)
       , prefix(needle.size())
   {
      for (size_t idx = 1; idx < this->needle.size(); idx++) {
         size_t border = prefix[idx - 1];
         while (border > 0 && this->needle[idx] != this->needle[border]) {
            border = prefix[border - 1];
         }
         if (this->needle[idx] == this->needle[border]) {
            border++;
         }
         prefix[idx] = border;
      }
   }

   void Reset()
   {
      state = 0;
      matched = needle.empty();
   }

   void Consume(std::string_view text)
   {
      if (matched) {
         return;
      }
      for (char c : text) {
         while (state > 0 && c != needle[state]) {
            state = prefix[state - 1];
         }
         if (c == needle[state]) {
            state++;
            if (state == needle.size()) {
               matched = true;
               return;
            }
         }
      }
   }

   bool Accepted() const
   {
      return matched;
   }

private:
   std::string needle;
   std::vector<size_t> prefix;
   size_t state = 0;
   bool matched = false;
};
// -------------------------------------------------------------------------------------
class TokenLikeEngine : public Engine {
public:
   enum class Mode {
      Contains,
      General
   };

   TokenLikeEngine(std::string_view pattern, Mode mode)
       : pattern(pattern)
       , mode(mode)
       , contains_matcher(mode == Mode::Contains ? pattern : std::string_view())
       , like_matcher(pattern)
   {
   }

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) final
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (MatchesText(block.GetRow(row_idx))) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const CompressedBlock& block, std::vector<uint32_t>& result) final
   {
      if (!block.IsTokenizerCodec()) {
         return 0;
      }

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (MatchesCompressed(block.GetRow(row_idx))) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

private:
   bool MatchesText(std::string_view text)
   {
      if (mode == Mode::Contains) {
         return text.find(pattern) != std::string_view::npos;
      }

      like_matcher.Reset();
      like_matcher.Consume(text);
      return like_matcher.Accepted();
   }

   bool MatchesCompressed(std::string_view compressed_text)
   {
      const size_t footer_size = sizeof(uint64_t);
      if (compressed_text.size() < footer_size) {
         return false;
      }

      const size_t token_bytes = compressed_text.size() - footer_size;
      if ((token_bytes % sizeof(uint16_t)) != 0) {
         return false;
      }

      auto consume_tokens = [&](auto& matcher) {
         matcher.Reset();
         const auto& dictionary = tokenizer_codec::TokenizerDictionary::Get();
         const char* token_data = compressed_text.data();
         for (size_t offset = 0; offset < token_bytes; offset += sizeof(uint16_t)) {
            uint16_t token_id = 0;
            std::memcpy(&token_id, token_data + offset, sizeof(token_id));
            matcher.Consume(dictionary.TokenText(token_id));
            if (matcher.Accepted() && mode == Mode::Contains) {
               return true;
            }
         }
         return matcher.Accepted();
      };

      if (mode == Mode::Contains) {
         return consume_tokens(contains_matcher);
      }
      return consume_tokens(like_matcher);
   }

   std::string pattern;
   Mode mode;
   TokenContainsMatcher contains_matcher;
   PercentLikeMatcher like_matcher;
};
// -------------------------------------------------------------------------------------
class TokenLikeEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (pattern.find('_') != std::string_view::npos) {
         return nullptr;
      }

      if (std::count(pattern.begin(), pattern.end(), '%') == 2 &&
          pattern.starts_with('%') &&
          pattern.ends_with('%') &&
          pattern.substr(1, pattern.size() - 2).find('%') == std::string_view::npos) {
         return std::make_unique<TokenLikeEngine>(pattern.substr(1, pattern.size() - 2),
                                                  TokenLikeEngine::Mode::Contains);
      }

      return std::make_unique<TokenLikeEngine>(pattern, TokenLikeEngine::Mode::General);
   }

   std::string GetName() final { return "token_like"; }
};
// -------------------------------------------------------------------------------------
