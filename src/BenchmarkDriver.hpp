#pragma once
// -------------------------------------------------------------------------------------
#include <bitset>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include "FsstWrapper.hpp"
#include "OnPairWrapper.hpp"
#include "Utility.hpp"
#include "fsst/fsst.h"
// -------------------------------------------------------------------------------------
constexpr uint32_t BLOCK_SIZE = 64 * 1024;
// -------------------------------------------------------------------------------------
struct FsstBlock {
   uint32_t row_count;
   std::vector<char> data;
   FsstDecoder decoder;
   std::bitset<256> used_chars;
   std::array<uint32_t, BLOCK_SIZE + 1> offsets;

   std::string_view GetRow(uint32_t row_idx) const
   {
      uint32_t start = offsets[row_idx];
      uint32_t end = offsets[row_idx + 1];
      return std::string_view(data.data() + start, end - start);
   }

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
struct OnPairBlock {
   uint32_t row_count = 0;
   uint32_t decompressed_bytes = 0;   // total raw bytes; lets bulk decoders size buffers
   OnPairColumnHandle column;
   std::array<uint32_t, BLOCK_SIZE + 1> offsets{};
};
// -------------------------------------------------------------------------------------
class Engine {
public:
   virtual ~Engine() = default;

   // Each overload defaults to throwing logic_error so engines opt in to the
   // block types they support. The driver catches and reports these.
   virtual uint32_t Scan(const RawBlock&, std::vector<uint32_t>&)
   {
      throw std::logic_error("Engine::Scan(RawBlock) not implemented");
   }
   virtual uint32_t Scan(const FsstBlock&, std::vector<uint32_t>&)
   {
      throw std::logic_error("Engine::Scan(FsstBlock) not implemented");
   }
   virtual uint32_t Scan(const OnPairBlock&, std::vector<uint32_t>&)
   {
      throw std::logic_error("Engine::Scan(OnPairBlock) not implemented");
   }
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
   void LoadBlocks(std::string_view file_path);
   void Run(std::string_view pattern);

private:
   std::vector<std::unique_ptr<EngineFactory>> engine_factories;
   std::vector<RawBlock> raw_blocks;
   std::vector<FsstBlock> fsst_blocks;
   std::vector<OnPairBlock> onpair_blocks;

   FsstBlock CreateFsstBlock(const RawBlock& raw_block) const;
   OnPairBlock CreateOnPairBlock(const RawBlock& raw_block) const;
};
// -------------------------------------------------------------------------------------
