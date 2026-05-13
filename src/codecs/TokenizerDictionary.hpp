#pragma once
// -------------------------------------------------------------------------------------
#include <cstdint>
#include <string_view>
#include <vector>
#include <onpair/core/dictionary_view.h>
// -------------------------------------------------------------------------------------
namespace tokenizer_codec {
// -------------------------------------------------------------------------------------
struct TokenRankRange {
   uint32_t begin = 1;
   uint32_t last = 0;

   bool Empty() const { return begin > last; }
   bool Contains(uint32_t rank) const { return rank >= begin && rank <= last; }
};
// -------------------------------------------------------------------------------------
class TokenizerDictionary {
public:
   static const TokenizerDictionary& Get();

   std::string_view TokenText(uint16_t token_id) const;
   std::string_view TokenTextByRank(uint32_t rank) const;
   uint32_t SortedRank(uint16_t token_id) const;
   uint16_t TokenIdByRank(uint32_t rank) const;
   TokenRankRange PrefixRange(const uint8_t* prefix, size_t prefix_len) const;
   onpair::DictionaryView OnPairDictionaryView() const { return onpair::DictionaryView(onpair_dictionary); }
   size_t TokenCount() const { return lengths.size(); }

private:
   TokenizerDictionary();

   std::vector<uint8_t> lengths;
   std::vector<uint64_t> tokens;
   std::vector<uint16_t> sorted_token_ids;
   std::vector<uint32_t> sorted_rank_by_token_id;
   onpair::Dictionary onpair_dictionary;
};
// -------------------------------------------------------------------------------------
} // namespace tokenizer_codec
// -------------------------------------------------------------------------------------
