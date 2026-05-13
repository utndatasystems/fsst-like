#include "BenchmarkDriver.hpp"
#include <chrono>
#include <fstream>
#include <stdexcept>
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
void BenchmarkDriver::AddEngine(unique_ptr<EngineFactory> engine_factory)
{
   engine_factories.push_back(move(engine_factory));
}
// -------------------------------------------------------------------------------------
void BenchmarkDriver::LoadBlocks(string_view file_path)
{
   raw_blocks.clear();
   fsst_blocks.clear();
   onpair_blocks.clear();

   // Open the file.
   string path_as_string(file_path);
   ifstream in(path_as_string);
   if (!in.is_open()) {
      throw runtime_error("Unable to open file '" + path_as_string + "'.");
   }

   // Read file line by line and insert into block.
   RawBlock block;
   block.offsets[0] = 0;
   string line;
   block.row_count = 0;
   while (getline(in, line)) {
      block.data.insert(block.data.end(), line.begin(), line.end());
      block.offsets[block.row_count + 1] = block.data.size();

      // If block is full -> append to raw_blocks and reset block.
      block.row_count++;
      if (block.row_count == BLOCK_SIZE) {
         raw_blocks.push_back(move(block));
         block.data = {};
         block.offsets[0] = 0;
         block.row_count = 0;
      }
   }

   // Do not forget the unfull block ;).
   if (block.row_count > 0) {
      raw_blocks.push_back(move(block));
   }

   // Compress all blocks (FSST).
   for (RawBlock& raw_block : raw_blocks) {
      fsst_blocks.push_back(CreateFsstBlock(raw_block));
   }

   // Compress all blocks (OnPair).
   onpair_blocks.reserve(raw_blocks.size());
   for (RawBlock& raw_block : raw_blocks) {
      onpair_blocks.push_back(CreateOnPairBlock(raw_block));
   }
}
// -------------------------------------------------------------------------------------
void BenchmarkDriver::Run(string_view pattern)
{
   vector<uint32_t> result(BLOCK_SIZE);
   for (auto& engine_factory : engine_factories) {
      // Create engine
      auto engine = engine_factory->Create(pattern);
      if (!engine) {
         std::cout << engine_factory->GetName() << " skipped (factory)" << std::endl;
         continue;
      }

      auto time_run = [&](auto& blocks, const char* label,
                          uint32_t& hits, long& ms) -> bool {
         try {
            auto t0 = std::chrono::high_resolution_clock::now();
            for (auto& b : blocks) {
               hits += engine->Scan(b, result);
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
            return true;
         } catch (const std::logic_error& e) {
            std::cout << engine_factory->GetName() << " " << label
                      << " skipped: " << e.what() << std::endl;
            return false;
         }
      };

      uint32_t raw_hits = 0, fsst_hits = 0, op_hits = 0;
      long raw_ms = -1, fsst_ms = -1, op_ms = -1;
      bool ok_raw  = time_run(raw_blocks,    "raw",    raw_hits,  raw_ms);
      bool ok_fsst = time_run(fsst_blocks,   "fsst",   fsst_hits, fsst_ms);
      bool ok_op   = time_run(onpair_blocks, "onpair", op_hits,   op_ms);

      std::cout << engine_factory->GetName()
                << ", raw="    << (ok_raw  ? std::to_string(raw_hits)  : "-")
                << " (" << raw_ms  << "ms)"
                << ", fsst="   << (ok_fsst ? std::to_string(fsst_hits) : "-")
                << " (" << fsst_ms << "ms)"
                << ", onpair=" << (ok_op   ? std::to_string(op_hits)   : "-")
                << " (" << op_ms   << "ms)" << endl;
   }
}
// -------------------------------------------------------------------------------------
FsstBlock BenchmarkDriver::CreateFsstBlock(const RawBlock& raw_block) const
{
   uint32_t row_count = raw_block.row_count;

   // Create row length and pointer array for fsst encoder.
   vector<const char*> ptrs(row_count);
   vector<uint64_t> lengths(row_count);
   for (uint32_t idx = 0; idx < row_count; idx++) {
      ptrs[idx] = raw_block.data.data() + raw_block.offsets[idx];
      lengths[idx] = raw_block.offsets[idx + 1] - raw_block.offsets[idx];
   }

   // Encode the data.
   FsstEncoder encoder;
   encoder.InitializeEncoderOnSample(row_count, (const char**)ptrs.data(), lengths.data());
   vector<char> compressed_buffer(7 + 2 * raw_block.data.size());
   vector<char*> compressed_ptrs(row_count);
   vector<uint64_t> compressed_lengths(row_count);
   uint32_t compressed_row_count = encoder.EncodeData(row_count,
                                                      (const char**)ptrs.data(), lengths.data(),          // in
                                                      compressed_buffer,                                  // out
                                                      compressed_ptrs.data(), compressed_lengths.data()); // out
   assert(compressed_row_count == row_count);

   // Store the compressed data in the block.
   FsstBlock fsst_block;
   fsst_block.row_count = row_count;
   uint32_t compressed_size = encoder.GetEncodedSize(row_count, compressed_ptrs.data(), compressed_lengths.data());
   fsst_block.data.resize(compressed_size);
   fsst_block.data.reserve(compressed_size + 128); // Reserve some space for easy SIMD.
   memcpy(fsst_block.data.data(), compressed_buffer.data(), compressed_size);
   for (uint32_t idx = 0; idx < row_count; idx++) {
      fsst_block.offsets[idx] = compressed_ptrs[idx] - compressed_buffer.data();
   }
   fsst_block.offsets[row_count] = compressed_size;

   // Create decoder.
   vector<char> decoder_buffer(encoder.GetRequiredDecoderSize());
   encoder.SerializeDecoder(decoder_buffer);
   fsst_block.decoder.DeserializeDecoder(decoder_buffer);

   // Find characters that occur in the encoded data not hidden behind symbols.
   uint8_t prev = 0;
   fsst_block.used_chars.reset();
   for (char c : fsst_block.data) {
      uint8_t current = static_cast<uint8_t>(c);
      if (prev == 255) {
         fsst_block.used_chars.set(current);
      }
      prev = current;
   }

   return fsst_block;
}
// -------------------------------------------------------------------------------------
OnPairBlock BenchmarkDriver::CreateOnPairBlock(const RawBlock& raw_block) const
{
   OnPairBlock blk;
   blk.row_count = raw_block.row_count;
   blk.decompressed_bytes = static_cast<uint32_t>(raw_block.data.size());
   blk.offsets = raw_block.offsets;

   onpair::encoding::TrainingConfig cfg;
   cfg.bits = 12;
   cfg.threshold = onpair::encoding::DynamicThreshold{1.0};
   cfg.seed = 42;

   blk.column.Build(raw_block.data.data(), raw_block.offsets.data(),
                    raw_block.row_count, cfg);
   return blk;
}
// -------------------------------------------------------------------------------------
