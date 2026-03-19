#pragma once
// -------------------------------------------------------------------------------------
#include "BenchmarkDriver.hpp"
// -------------------------------------------------------------------------------------
class PassThroughEngine : public Engine {
public:
   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) final
   {
      return ScanBlock(block, result);
   }

   uint32_t Scan(const CompressedBlock& block, std::vector<uint32_t>& result) final
   {
      return ScanBlock(block, result);
   }

private:
   template <class BlockType>
   uint32_t ScanBlock(const BlockType& block, std::vector<uint32_t>& result)
   {
      uint32_t count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         std::string_view row = block.GetRow(row_idx);
         checksum += static_cast<uint64_t>(row.size());
         result[count++] = row_idx;
      }
      return count;
   }

   volatile uint64_t checksum = 0;
};
// -------------------------------------------------------------------------------------
class PassThroughEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view /*pattern*/) final
   {
      return std::make_unique<PassThroughEngine>();
   }

   std::string GetName() final { return "pass_through"; }
};
// -------------------------------------------------------------------------------------
