#include "codecs/OnPairPlusCodec.hpp"
#include "compressor/onpair_advanced/OnPairAdvancedCompressor.hpp"

#include <cstddef>

namespace onpairplus_codec {

namespace {

sgtt::compressor::OnPairAdvancedCompressor<sgtt::compressor::onpair::MaxSymbolLength::SIXTEEN>& GetOnPairPlusCompressor()
{
   static sgtt::compressor::OnPairConfig config{};
   static sgtt::compressor::OnPairAdvancedCompressor<sgtt::compressor::onpair::MaxSymbolLength::SIXTEEN> compressor(config);
   return compressor;
}

} // namespace

void Compress(std::span<const char> input, std::vector<char>& output)
{
   auto& compressor = GetOnPairPlusCompressor();

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
   auto& compressor = GetOnPairPlusCompressor();

   std::span<std::byte> in_bytes(reinterpret_cast<std::byte*>(const_cast<char*>(input.data())), input.size());
   const size_t required_size = compressor.getUncompressedSize(in_bytes);
   std::vector<std::byte> decompressed(required_size);
   compressor.decompress(in_bytes, decompressed);

   output.resize(decompressed.size());
   for (size_t idx = 0; idx < decompressed.size(); ++idx) {
      output[idx] = static_cast<char>(decompressed[idx]);
   }
}

} // namespace onpairplus_codec

