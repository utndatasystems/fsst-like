#include "decoders/OnPairDecoder.hpp"
#include "codecs/OnPairCodec.hpp"

#include <cstring>
#include <vector>

uint32_t OnPairDecoder::Decode(std::span<const char> input, std::span<char> output) const
{
   std::vector<char> decompressed;
   onpair_codec::Decompress(input, decompressed);
   if (decompressed.size() > output.size()) {
      return decompressed.size();
   }

   std::memcpy(output.data(), decompressed.data(), decompressed.size());
   return static_cast<uint32_t>(decompressed.size());
}

uint32_t OnPairDecoder::GetIdealBufferSize(uint32_t compressed_size) const
{
   return compressed_size * 16 + 64;
}
