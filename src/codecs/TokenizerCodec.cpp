#include "codecs/TokenizerCodec.hpp"
#include "compressor/TokenizerCompressor.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace tokenizer_codec {

namespace {

sgtt::compressor::TokenizerCompressor& GetTokenizerCompressor()
{
   static sgtt::compressor::TokenizerCompressor compressor;
   return compressor;
}

} // namespace

void Compress(std::span<const char> input, std::vector<char>& output)
{
   auto& compressor = GetTokenizerCompressor();

   std::span<std::byte> in_bytes(reinterpret_cast<std::byte*>(const_cast<char*>(input.data())), input.size());
   std::vector<std::byte> compressed;
   compressor.compress(in_bytes, compressed);

   output.resize(compressed.size());
   for (size_t idx = 0; idx < compressed.size(); ++idx) {
      output[idx] = static_cast<char>(compressed[idx]);
   }
}

void Decompress(std::span<const char> input, std::vector<char>& output)
{
   auto& compressor = GetTokenizerCompressor();

   std::span<std::byte> in_bytes(reinterpret_cast<std::byte*>(const_cast<char*>(input.data())), input.size());
   std::vector<std::byte> decompressed;
   compressor.decompress(in_bytes, decompressed);

   output.resize(decompressed.size());
   for (size_t idx = 0; idx < decompressed.size(); ++idx) {
      output[idx] = static_cast<char>(decompressed[idx]);
   }
}

} // namespace tokenizer_codec
