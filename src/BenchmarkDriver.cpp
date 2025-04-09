#include "BenchmarkDriver.hpp"
#include <chrono>
#include <fstream>
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

   // Compress all blocks.
   for (RawBlock& raw_block : raw_blocks) {
      fsst_blocks.push_back(CreateFsstBlock(raw_block));
   }
}
// -------------------------------------------------------------------------------------
void BenchmarkDriver::RunRaw(string_view pattern)
{
   vector<uint32_t> result(BLOCK_SIZE);
   for (auto& engine_factory : engine_factories) {
      std::cout << "Running raw scan with " << engine_factory->GetName() << std::endl;
      uint32_t row_count = 0;
      auto engine = engine_factory->Create(pattern);
      auto begin = std::chrono::high_resolution_clock::now();
      for (auto& block : raw_blocks) {
         row_count += engine->Scan(block, result);
      }
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
      cout << "Found " << row_count << " in " << duration << "ms" << endl;
   }
}
// -------------------------------------------------------------------------------------
void BenchmarkDriver::RunCompressed(string_view pattern)
{
   vector<uint32_t> result(BLOCK_SIZE);
   for (auto& engine_factory : engine_factories) {
      std::cout << "Running raw scan with " << engine_factory->GetName() << std::endl;
      uint32_t row_count = 0;
      auto engine = engine_factory->Create(pattern);
      auto begin = std::chrono::high_resolution_clock::now();
      for (auto& block : fsst_blocks) {
         row_count += engine->Scan(block, result);
      }
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
      cout << "Found " << row_count << " in " << duration << "ms" << endl;
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
