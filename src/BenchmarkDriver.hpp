#pragma once
// -------------------------------------------------------------------------------------
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "FsstWrapper.hpp"
#include "Utility.hpp"
#include "fsst/fsst.h"
// -------------------------------------------------------------------------------------
constexpr uint32_t BLOCK_SIZE = 64 * 1024;
// -------------------------------------------------------------------------------------
struct FsstBlock {
   uint32_t row_count;
   std::vector<char> data;
   std::array<uint32_t, BLOCK_SIZE + 1> offsets;
   FsstDecoder decoder;
};
// -------------------------------------------------------------------------------------
struct RawBlock {
   uint32_t row_count;
   std::vector<char> data;
   std::array<uint32_t, BLOCK_SIZE + 1> offsets;
};
// -------------------------------------------------------------------------------------
class Engine {
public:
   // Called once for each block.
   virtual uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result) = 0;
   virtual uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) = 0;
};
// -------------------------------------------------------------------------------------
class EngineFactory {
public:
   // Called once per table, at the begining of the scan operation.
   virtual std::unique_ptr<Engine> Create(std::string_view pattern) = 0;
   virtual std::string GetName() = 0;
};
// -------------------------------------------------------------------------------------
class BenchmarkDriver {
public:
   void AddEngine(std::unique_ptr<EngineFactory> engine_factory);
   void LoadBlocks(std::string_view file_path);
   void RunRaw(std::string_view pattern);
   void RunCompressed(std::string_view pattern);

private:
   std::vector<std::unique_ptr<EngineFactory>> engine_factories;
   std::vector<RawBlock> raw_blocks;
   std::vector<FsstBlock> fsst_blocks;
};
// -------------------------------------------------------------------------------------
