#include "decoders/TokenizerDecoder.hpp"
#include "codecs/TokenizerCodec.hpp"

#include <cstring>

uint32_t TokenizerDecoder::Decode(std::span<const char> input, std::span<char> output) const
{
   std::vector<char> decompressed;
   tokenizer_codec::Decompress(input, decompressed);
   if (decompressed.size() > output.size()) {
      return decompressed.size();
   }

   std::memcpy(output.data(), decompressed.data(), decompressed.size());
   return static_cast<uint32_t>(decompressed.size());
}

uint32_t TokenizerDecoder::GetIdealBufferSize(uint32_t compressed_size) const
{
   return compressed_size * 8 + 32;
}
