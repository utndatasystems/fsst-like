#pragma once

#include <span>
#include <vector>

namespace onpairplus_codec {

void Compress(std::span<const char> input, std::vector<char>& output);
void Decompress(std::span<const char> input, std::vector<char>& output);

} // namespace onpairplus_codec

