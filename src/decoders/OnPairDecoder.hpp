#pragma once

#include <cstdint>
#include <span>

class OnPairDecoder {
public:
   uint32_t Decode(std::span<const char> input, std::span<char> output) const;
   uint32_t GetIdealBufferSize(uint32_t compressed_size) const;
};

