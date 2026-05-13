#include "codecs/TokenizerDictionary.hpp"

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
} // namespace tokenizer_codec
// -------------------------------------------------------------------------------------
