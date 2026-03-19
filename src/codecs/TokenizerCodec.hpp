#pragma once

#include <span>
#include <vector>

namespace tokenizer_codec {

void Compress(std::span<const char> input, std::vector<char>& output);
void Decompress(std::span<const char> input, std::vector<char>& output);

} // namespace tokenizer_codec

