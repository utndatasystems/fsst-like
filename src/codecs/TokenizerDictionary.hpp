#pragma once
// -------------------------------------------------------------------------------------
#include <cstdint>
#include <string_view>
#include <vector>
// -------------------------------------------------------------------------------------
namespace tokenizer_codec {
// -------------------------------------------------------------------------------------
class TokenizerDictionary {
public:
   static const TokenizerDictionary& Get();

   std::string_view TokenText(uint16_t token_id) const;
   size_t TokenCount() const { return lengths.size(); }

private:
   TokenizerDictionary();

   std::vector<uint8_t> lengths;
   std::vector<uint64_t> tokens;
};
// -------------------------------------------------------------------------------------
} // namespace tokenizer_codec
// -------------------------------------------------------------------------------------
