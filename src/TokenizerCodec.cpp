#include "TokenizerCodec.hpp"
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>

using namespace std;

namespace tokenizer_codec {

namespace {

struct TrieNode {
   int token_id = -1;
   unordered_map<uint8_t, uint32_t> next;
};

struct Dictionary {
   vector<vector<uint8_t>> tokens;
   vector<TrieNode> trie;
};

uint8_t HexNibble(char c)
{
   if (c >= '0' && c <= '9') return static_cast<uint8_t>(c - '0');
   if (c >= 'a' && c <= 'f') return static_cast<uint8_t>(10 + (c - 'a'));
   if (c >= 'A' && c <= 'F') return static_cast<uint8_t>(10 + (c - 'A'));
   throw runtime_error("Invalid hex character in token CSV.");
}

vector<uint8_t> ParseHexToken(string_view hex)
{
   if ((hex.size() % 2) != 0) {
      throw runtime_error("Invalid hex token length in token CSV.");
   }

   vector<uint8_t> bytes;
   bytes.reserve(hex.size() / 2);
   for (size_t idx = 0; idx < hex.size(); idx += 2) {
      uint8_t high = HexNibble(hex[idx]);
      uint8_t low = HexNibble(hex[idx + 1]);
      bytes.push_back(static_cast<uint8_t>((high << 4) | low));
   }
   return bytes;
}

const Dictionary& GetDictionary()
{
   static once_flag once;
   static Dictionary dict;
   call_once(once, [&] {
      ifstream in("third-party/token-vldb2026/tokens/openai_100k_patched.csv");
      if (!in.is_open()) {
         throw runtime_error("Unable to open token dictionary CSV.");
      }

      dict.trie.push_back({});

      string line;
      getline(in, line); // header
      while (getline(in, line)) {
         size_t pos1 = line.find(',');
         if (pos1 == string::npos) continue;
         size_t pos2 = line.find(',', pos1 + 1);
         if (pos2 == string::npos) continue;

         string token_hex = line.substr(0, pos1);
         string rank_str = line.substr(pos1 + 1, pos2 - (pos1 + 1));

         uint32_t token_id = stoul(rank_str);
         vector<uint8_t> token_bytes = ParseHexToken(token_hex);

         if (dict.tokens.size() <= token_id) {
            dict.tokens.resize(token_id + 1);
         }
         dict.tokens[token_id] = token_bytes;

         uint32_t node_idx = 0;
         for (uint8_t byte : token_bytes) {
            auto [it, inserted] = dict.trie[node_idx].next.emplace(byte, 0);
            if (inserted) {
               it->second = dict.trie.size();
               dict.trie.push_back({});
            }
            node_idx = it->second;
         }
         dict.trie[node_idx].token_id = token_id;
      }
   });

   return dict;
}

} // namespace

void Compress(std::span<const char> input, std::vector<char>& output)
{
   const Dictionary& dict = GetDictionary();

   output.clear();
   output.reserve(input.size() * sizeof(uint16_t) + sizeof(uint64_t));

   const uint8_t* bytes = reinterpret_cast<const uint8_t*>(input.data());
   size_t idx = 0;
   while (idx < input.size()) {
      int best_token = -1;
      size_t best_len = 0;

      uint32_t node_idx = 0;
      for (size_t jdx = idx; jdx < input.size(); jdx++) {
         uint8_t byte = bytes[jdx];
         auto it = dict.trie[node_idx].next.find(byte);
         if (it == dict.trie[node_idx].next.end()) break;
         node_idx = it->second;
         if (dict.trie[node_idx].token_id >= 0) {
            best_token = dict.trie[node_idx].token_id;
            best_len = jdx - idx + 1;
         }
      }

      if (best_token < 0) {
         best_token = bytes[idx];
         best_len = 1;
      }

      uint16_t token_id = static_cast<uint16_t>(best_token);
      output.push_back(static_cast<char>(token_id & 0xFF));
      output.push_back(static_cast<char>((token_id >> 8) & 0xFF));
      idx += best_len;
   }

   uint64_t uncompressed_size = input.size();
   size_t old_size = output.size();
   output.resize(old_size + sizeof(uint64_t));
   memcpy(output.data() + old_size, &uncompressed_size, sizeof(uint64_t));
}

void Decompress(std::span<const char> input, std::vector<char>& output)
{
   const Dictionary& dict = GetDictionary();

   if (input.size() < sizeof(uint64_t)) {
      throw runtime_error("Tokenizer input too small.");
   }

   size_t token_bytes = input.size() - sizeof(uint64_t);
   if ((token_bytes % sizeof(uint16_t)) != 0) {
      throw runtime_error("Tokenizer token stream has odd byte size.");
   }

   uint64_t uncompressed_size = 0;
   memcpy(&uncompressed_size, input.data() + token_bytes, sizeof(uint64_t));

   output.clear();
   output.reserve(uncompressed_size);

   const uint8_t* bytes = reinterpret_cast<const uint8_t*>(input.data());
   for (size_t idx = 0; idx < token_bytes; idx += 2) {
      uint16_t token_id = static_cast<uint16_t>(bytes[idx] | (static_cast<uint16_t>(bytes[idx + 1]) << 8));
      if (token_id >= dict.tokens.size()) {
         throw runtime_error("Tokenizer token id out of range.");
      }
      const auto& token = dict.tokens[token_id];
      output.insert(output.end(), token.begin(), token.end());
   }
}

} // namespace tokenizer_codec

