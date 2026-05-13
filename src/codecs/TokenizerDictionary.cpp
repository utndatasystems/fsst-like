#include "codecs/TokenizerDictionary.hpp"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
// -------------------------------------------------------------------------------------
namespace tokenizer_codec {
// -------------------------------------------------------------------------------------
namespace {
// -------------------------------------------------------------------------------------
uint64_t LoadU64Le(const uint8_t* data)
{
   uint64_t value = 0;
   for (size_t idx = 0; idx < sizeof(uint64_t); idx++) {
      value |= static_cast<uint64_t>(data[idx]) << (idx * 8);
   }
   return value;
}
// -------------------------------------------------------------------------------------
std::vector<std::filesystem::path> BinCandidates(const char* file_name)
{
   namespace fs = std::filesystem;
   std::vector<fs::path> paths;
#if defined(SGTT_TOKENIZER_BIN_DIR)
   paths.emplace_back(fs::path(SGTT_TOKENIZER_BIN_DIR) / file_name);
#endif
   paths.emplace_back(fs::path("src/tokens/openai100kpatched") / file_name);
   paths.emplace_back(fs::path("../src/tokens/openai100kpatched") / file_name);
   paths.emplace_back(fs::path("../../src/tokens/openai100kpatched") / file_name);
   paths.emplace_back(fs::path("third-party/token-vldb2026/src/tokens/openai100kpatched") / file_name);
   paths.emplace_back(fs::path("../third-party/token-vldb2026/src/tokens/openai100kpatched") / file_name);
   return paths;
}
// -------------------------------------------------------------------------------------
std::vector<uint8_t> ReadFirstExistingFile(const char* file_name)
{
   for (const auto& path : BinCandidates(file_name)) {
      std::ifstream file(path, std::ios::binary);
      if (!file.is_open()) {
         continue;
      }

      file.seekg(0, std::ios::end);
      const size_t size = static_cast<size_t>(file.tellg());
      file.seekg(0, std::ios::beg);

      std::vector<uint8_t> bytes(size);
      file.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
      if (!file) {
         throw std::runtime_error("TokenizerDictionary: failed reading " + path.string());
      }
      return bytes;
   }
   throw std::runtime_error(std::string("TokenizerDictionary: unable to open ") + file_name);
}
// -------------------------------------------------------------------------------------
} // namespace
// -------------------------------------------------------------------------------------
const TokenizerDictionary& TokenizerDictionary::Get()
{
   static const TokenizerDictionary dictionary;
   return dictionary;
}
// -------------------------------------------------------------------------------------
TokenizerDictionary::TokenizerDictionary()
{
   std::vector<uint8_t> runtime_lengths = ReadFirstExistingFile("openai100kpatched_lengths.bin");
   std::vector<uint8_t> runtime_dict = ReadFirstExistingFile("openai100kpatched_dict.bin");

   if (runtime_lengths.empty() || runtime_dict.size() % runtime_lengths.size() != 0) {
      throw std::runtime_error("TokenizerDictionary: invalid runtime dictionary files");
   }

   const size_t token_count = runtime_lengths.size();
   const size_t bytes_per_token = runtime_dict.size() / token_count;
   if (bytes_per_token < sizeof(uint64_t) || token_count > 65536) {
      throw std::runtime_error("TokenizerDictionary: unsupported runtime dictionary layout");
   }

   lengths.resize(token_count);
   tokens.resize(token_count);
   for (size_t token_id = 0; token_id < token_count; token_id++) {
      const uint8_t length = runtime_lengths[token_id];
      if (length == 0 || length > sizeof(uint64_t)) {
         throw std::runtime_error("TokenizerDictionary: invalid token length");
      }
      lengths[token_id] = length;
      tokens[token_id] = LoadU64Le(runtime_dict.data() + token_id * bytes_per_token);
   }

   sorted_token_ids.resize(token_count);
   sorted_rank_by_token_id.resize(token_count);
   for (size_t token_id = 0; token_id < token_count; token_id++) {
      sorted_token_ids[token_id] = static_cast<uint16_t>(token_id);
   }

   std::sort(sorted_token_ids.begin(), sorted_token_ids.end(), [&](uint16_t lhs, uint16_t rhs) {
      const auto left = TokenText(lhs);
      const auto right = TokenText(rhs);
      const int cmp = std::memcmp(left.data(), right.data(), std::min(left.size(), right.size()));
      if (cmp != 0) {
         return cmp < 0;
      }
      if (left.size() != right.size()) {
         return left.size() < right.size();
      }
      return lhs < rhs;
   });

   for (uint32_t rank = 0; rank < sorted_token_ids.size(); rank++) {
      sorted_rank_by_token_id[sorted_token_ids[rank]] = rank;
   }

   onpair_dictionary.offsets.reserve(token_count + 1);
   onpair_dictionary.offsets.push_back(0);
   for (uint16_t token_id : sorted_token_ids) {
      const auto text = TokenText(token_id);
      onpair_dictionary.bytes.insert(onpair_dictionary.bytes.end(), text.begin(), text.end());
      onpair_dictionary.offsets.push_back(static_cast<uint32_t>(onpair_dictionary.bytes.size()));
   }
   onpair_dictionary.pad_for_decoder();
}
// -------------------------------------------------------------------------------------
std::string_view TokenizerDictionary::TokenText(uint16_t token_id) const
{
   if (token_id >= lengths.size()) {
      throw std::runtime_error("TokenizerDictionary: token id out of range");
   }
   return std::string_view(reinterpret_cast<const char*>(&tokens[token_id]), lengths[token_id]);
}
// -------------------------------------------------------------------------------------
std::string_view TokenizerDictionary::TokenTextByRank(uint32_t rank) const
{
   if (rank >= sorted_token_ids.size()) {
      throw std::runtime_error("TokenizerDictionary: sorted rank out of range");
   }
   return TokenText(sorted_token_ids[rank]);
}
// -------------------------------------------------------------------------------------
uint32_t TokenizerDictionary::SortedRank(uint16_t token_id) const
{
   if (token_id >= sorted_rank_by_token_id.size()) {
      throw std::runtime_error("TokenizerDictionary: token id out of range");
   }
   return sorted_rank_by_token_id[token_id];
}
// -------------------------------------------------------------------------------------
uint16_t TokenizerDictionary::TokenIdByRank(uint32_t rank) const
{
   if (rank >= sorted_token_ids.size()) {
      throw std::runtime_error("TokenizerDictionary: sorted rank out of range");
   }
   return sorted_token_ids[rank];
}
// -------------------------------------------------------------------------------------
TokenRankRange TokenizerDictionary::PrefixRange(const uint8_t* prefix, size_t prefix_len) const
{
   if (prefix_len > sizeof(uint64_t)) {
      return {};
   }

   auto lower_bound = [&](const uint8_t* target, size_t target_len, uint32_t start) {
      uint32_t lo = start;
      uint32_t hi = static_cast<uint32_t>(sorted_token_ids.size());
      while (lo < hi) {
         const uint32_t mid = lo + ((hi - lo) >> 1);
         const auto token = TokenTextByRank(mid);
         const size_t cmp_len = std::min(token.size(), target_len);
         const int cmp = std::memcmp(token.data(), target, cmp_len);
         if (cmp < 0 || (cmp == 0 && token.size() < target_len)) {
            lo = mid + 1;
         } else {
            hi = mid;
         }
      }
      return lo;
   };

   const uint32_t lo = lower_bound(prefix, prefix_len, 0);

   uint8_t upper[sizeof(uint64_t)];
   size_t upper_len = prefix_len;
   bool overflow = true;
   while (upper_len > 0) {
      if (prefix[upper_len - 1] < 0xFF) {
         std::memcpy(upper, prefix, upper_len);
         upper[upper_len - 1]++;
         overflow = false;
         break;
      }
      upper_len--;
   }

   const uint32_t hi = overflow ? static_cast<uint32_t>(sorted_token_ids.size()) : lower_bound(upper, upper_len, lo);
   if (lo >= hi) {
      return {};
   }
   return TokenRankRange{lo, hi - 1};
}
// -------------------------------------------------------------------------------------
} // namespace tokenizer_codec
// -------------------------------------------------------------------------------------
