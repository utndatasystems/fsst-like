#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>
#include "codecs/TokenizerCodec.hpp"
#include "codecs/TokenizerDictionary.hpp"
// -------------------------------------------------------------------------------------
namespace token_like {
// -------------------------------------------------------------------------------------
inline std::vector<uint16_t> TokenizePattern(std::string_view pattern)
{
   if (pattern.empty()) {
      return {};
   }

   std::vector<char> compressed;
   tokenizer_codec::Compress(std::span<const char>(pattern.data(), pattern.size()), compressed);

   constexpr size_t footer_size = sizeof(uint64_t);
   if (compressed.size() < footer_size) {
      return {};
   }

   const size_t token_bytes = compressed.size() - footer_size;
   std::vector<uint16_t> tokens(token_bytes / sizeof(uint16_t));
   std::memcpy(tokens.data(), compressed.data(), token_bytes);
   return tokens;
}
// -------------------------------------------------------------------------------------
class TokenCursor {
public:
   explicit TokenCursor(std::string_view compressed_text)
       : data(compressed_text.data())
   {
      constexpr size_t footer_size = sizeof(uint64_t);
      if (compressed_text.size() < footer_size) {
         return;
      }
      token_bytes = compressed_text.size() - footer_size;
      if ((token_bytes % sizeof(uint16_t)) != 0) {
         token_bytes = 0;
      }
   }

   bool HasMore() const { return offset < token_bytes; }

   uint16_t Next()
   {
      uint16_t token = 0;
      std::memcpy(&token, data + offset, sizeof(token));
      offset += sizeof(token);
      return token;
   }

private:
   const char* data = nullptr;
   size_t token_bytes = 0;
   size_t offset = 0;
};
// -------------------------------------------------------------------------------------
template <class Automaton>
bool Drive(Automaton& automaton, std::string_view compressed_text)
{
   TokenCursor cursor(compressed_text);
   automaton.Reset();
   while (cursor.HasMore()) {
      automaton.Step(cursor.Next());
      if (automaton.IsDead()) {
         break;
      }
   }
   return automaton.IsAccepted();
}
// -------------------------------------------------------------------------------------
class KmpAutomaton {
public:
   using State = uint32_t;

   explicit KmpAutomaton(std::string_view pattern)
       : pattern(pattern)
       , match_state(static_cast<State>(pattern.size()))
   {
      Build();
   }

   void Step(uint16_t token) noexcept
   {
      if (IsDead()) {
         return;
      }
      state = transitions[(static_cast<size_t>(state) * token_count) + token];
   }

   bool IsAccepted() const noexcept { return state == match_state; }
   bool IsDead() const noexcept { return state == match_state; }
   void Reset() noexcept { state = 0; }

private:
   State StepBytes(State current, std::string_view bytes) const
   {
      for (char c : bytes) {
         if (current == match_state) {
            return match_state;
         }
         while (current > 0 && pattern[current] != c) {
            current = failure[current - 1];
         }
         if (pattern[current] == c) {
            current++;
         }
      }
      return current;
   }

   void Build()
   {
      const auto& dictionary = tokenizer_codec::TokenizerDictionary::Get();
      token_count = dictionary.TokenCount();
      transitions.assign((static_cast<size_t>(match_state) + 1) * token_count, 0);

      if (pattern.empty()) {
         return;
      }

      failure.assign(pattern.size(), 0);
      for (State idx = 1, border = 0; idx < pattern.size();) {
         if (pattern[idx] == pattern[border]) {
            failure[idx++] = ++border;
         } else if (border > 0) {
            border = failure[border - 1];
         } else {
            failure[idx++] = 0;
         }
      }

      for (State entry_state = 0; entry_state <= match_state; entry_state++) {
         for (uint32_t token = 0; token < token_count; token++) {
            State target = match_state;
            if (entry_state != match_state) {
               target = StepBytes(entry_state, dictionary.TokenText(static_cast<uint16_t>(token)));
            }
            transitions[(static_cast<size_t>(entry_state) * token_count) + token] = target;
         }
      }
   }

   std::string pattern;
   State match_state = 0;
   State state = 0;
   size_t token_count = 0;
   std::vector<State> failure;
   std::vector<State> transitions;
};
// -------------------------------------------------------------------------------------
class PrefixAutomaton {
public:
   explicit PrefixAutomaton(std::string_view prefix)
       : query_tokens(TokenizePattern(prefix))
   {
      const auto& dictionary = tokenizer_codec::TokenizerDictionary::Get();
      intervals.resize(query_tokens.size());

      size_t current_pos = 0;
      const auto* prefix_bytes = reinterpret_cast<const uint8_t*>(prefix.data());
      for (size_t idx = 0; idx < query_tokens.size(); idx++) {
         intervals[idx] = dictionary.PrefixRange(prefix_bytes + current_pos, prefix.size() - current_pos);
         current_pos += dictionary.TokenText(query_tokens[idx]).size();
      }
   }

   void Step(uint16_t token) noexcept
   {
      if (IsDead()) {
         return;
      }

      if (token != query_tokens[pos]) {
         const auto& dictionary = tokenizer_codec::TokenizerDictionary::Get();
         status = intervals[pos].Contains(dictionary.SortedRank(token)) ? Status::Accepted : Status::Rejected;
         return;
      }

      pos++;
      if (pos == query_tokens.size()) {
         status = Status::Accepted;
      }
   }

   bool IsAccepted() const noexcept { return status == Status::Accepted; }
   bool IsDead() const noexcept { return status != Status::Matching; }

   void Reset() noexcept
   {
      pos = 0;
      status = query_tokens.empty() ? Status::Accepted : Status::Matching;
   }

private:
   enum class Status : uint8_t {
      Matching,
      Accepted,
      Rejected
   };

   std::vector<uint16_t> query_tokens;
   std::vector<tokenizer_codec::TokenRankRange> intervals;
   size_t pos = 0;
   Status status = Status::Matching;
};
// -------------------------------------------------------------------------------------
} // namespace token_like
// -------------------------------------------------------------------------------------
