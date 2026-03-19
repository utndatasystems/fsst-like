#include "decoders/OnPairPlusDecoder.hpp"
#include "codecs/OnPairPlusCodec.hpp"

#include <cstring>
#include <vector>

uint32_t OnPairPlusDecoder::Decode(std::span<const char> input, std::span<char> output) const
{
   std::vector<char> decompressed;
   onpairplus_codec::Decompress(input, decompressed);
   if (decompressed.size() > output.size()) {
      return decompressed.size();
   }

   std::memcpy(output.data(), decompressed.data(), decompressed.size());
   return static_cast<uint32_t>(decompressed.size());
}

uint32_t OnPairPlusDecoder::GetIdealBufferSize(uint32_t compressed_size) const
{
   return compressed_size * 16 + 64;
}

