#pragma once
// -------------------------------------------------------------------------------------
#include <cstdint>
#include <cstring>
#include <string_view>
#include <onpair/search/automata/kmp_automaton.h>
#include <onpair/search/automata/prefix_automaton.h>
#include "codecs/TokenizerDictionary.hpp"
// -------------------------------------------------------------------------------------
namespace token_like {
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
class OnPairKmpAdapter {
public:
   explicit OnPairKmpAdapter(std::string_view pattern)
       : automaton(pattern, tokenizer_codec::TokenizerDictionary::Get().OnPairDictionaryView())
   {
   }

   void Step(uint16_t token) noexcept
   {
      automaton.step(static_cast<onpair::Token>(tokenizer_codec::TokenizerDictionary::Get().SortedRank(token)));
   }

   bool IsAccepted() const noexcept { return automaton.is_accepted(); }
   bool IsDead() const noexcept { return automaton.is_dead(); }
   void Reset() noexcept { automaton.reset(); }

private:
   onpair::search::KmpAutomaton automaton;
};
// -------------------------------------------------------------------------------------
class OnPairPrefixAdapter {
public:
   explicit OnPairPrefixAdapter(std::string_view prefix)
       : automaton(prefix, tokenizer_codec::TokenizerDictionary::Get().OnPairDictionaryView())
   {
   }

   void Step(uint16_t token) noexcept
   {
      automaton.step(static_cast<onpair::Token>(tokenizer_codec::TokenizerDictionary::Get().SortedRank(token)));
   }

   bool IsAccepted() const noexcept { return automaton.is_accepted(); }
   bool IsDead() const noexcept { return automaton.is_dead(); }
   void Reset() noexcept { automaton.reset(); }

private:
   onpair::search::PrefixAutomaton automaton;
};
} // namespace token_like
// -------------------------------------------------------------------------------------
