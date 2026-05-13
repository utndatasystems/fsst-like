#pragma once
// -------------------------------------------------------------------------------------
#include <bitset>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <onpair/api.h>
#include "decoders/FsstWrapper.hpp"
#include "decoders/TokenizerDecoder.hpp"
#include "Utility.hpp"
#include "fsst/fsst.h"
// -------------------------------------------------------------------------------------
constexpr uint32_t BLOCK_SIZE = 64 * 1024;
enum class CompressionCodec {
   Fsst,
   Tokenizer,
   OnPair
};
// -------------------------------------------------------------------------------------
struct CompressedBlock {
   uint32_t row_count;
   std::vector<char> data;
   CompressionCodec codec = CompressionCodec::Fsst;
   FsstDecoder fsst_decoder;
   TokenizerDecoder tokenizer_decoder;
   onpair::OnPairColumn onpair_column;
   uint32_t max_uncompressed_row_size = 0;
   std::bitset<256> used_chars;
   std::array<uint32_t, BLOCK_SIZE + 1> offsets;

   std::string_view GetRow(uint32_t row_idx) const
   {
      uint32_t start = offsets[row_idx];
      uint32_t end = offsets[row_idx + 1];
      return std::string_view(data.data() + start, end - start);
   }


   uint32_t Decode(std::string_view input, std::span<char> output) const
   {
      std::span<const char> encoded(input.data(), input.size());
      if (codec == CompressionCodec::Fsst) {
         return fsst_decoder.Decode(encoded, output);
      }
      if (codec == CompressionCodec::Tokenizer) {
         return tokenizer_decoder.Decode(encoded, output);
      }
      return 0;
   }

   uint32_t GetIdealBufferSize(uint32_t compressed_size) const
   {
      if (codec == CompressionCodec::Fsst) {
         return fsst_decoder.GetIdealBufferSize(compressed_size);
      }
      if (codec == CompressionCodec::Tokenizer) {
         return tokenizer_decoder.GetIdealBufferSize(compressed_size);
      }
      return max_uncompressed_row_size + onpair::DECOMPRESS_BUFFER_PADDING;
   }

   bool IsTokenizerCodec() const { return codec == CompressionCodec::Tokenizer; }
   bool IsFsstCodec() const { return codec == CompressionCodec::Fsst; }
   bool IsOnPairCodec() const { return codec == CompressionCodec::OnPair; }

   void PrintUsedChars(std::ostream& os) const
   {
      os << "used_chars: {" << std::endl;
      for (uint32_t idx = 0; idx < 256; idx++) {
         if (used_chars[idx]) {
            os << ((char)idx) << ", ";
         }
      }
      os << "}" << std::endl;
   }
};
// -------------------------------------------------------------------------------------
struct RawBlock {
   uint32_t row_count;
   std::vector<char> data;
   std::array<uint32_t, BLOCK_SIZE + 1> offsets;

   std::string_view GetRow(uint32_t row_idx) const
   {
      uint32_t start = offsets[row_idx];
      uint32_t end = offsets[row_idx + 1];
      return std::string_view(data.data() + start, end - start);
   }
};
// -------------------------------------------------------------------------------------
class Engine {
public:
   virtual ~Engine() = default;

   // Called once for each block.
   virtual uint32_t Scan(const CompressedBlock& block, std::vector<uint32_t>& result) = 0;
   virtual uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) = 0;
};
// -------------------------------------------------------------------------------------
class EngineFactory {
public:
   virtual ~EngineFactory() = default;

   // Called once per table, at the begining of the scan operation.
   virtual std::unique_ptr<Engine> Create(std::string_view pattern) = 0;
   virtual std::string GetName() = 0;
};
// -------------------------------------------------------------------------------------
class BenchmarkDriver {
public:
   void AddEngine(std::unique_ptr<EngineFactory> engine_factory);
   void LoadBlocks(std::string_view file_path, CompressionCodec codec = CompressionCodec::Fsst);
   void Run(std::string_view pattern);

private:
   std::vector<std::unique_ptr<EngineFactory>> engine_factories;
   std::vector<RawBlock> raw_blocks;
   std::vector<CompressedBlock> compressed_blocks;
   CompressionCodec codec = CompressionCodec::Fsst;

   CompressedBlock CreateFsstBlock(const RawBlock& raw_block) const;
   CompressedBlock CreateTokenizerBlock(const RawBlock& raw_block) const;
   CompressedBlock CreateOnPairBlock(const RawBlock& raw_block) const;
};
// -------------------------------------------------------------------------------------
