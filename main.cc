// this software is distributed under the MIT License (http://www.opensource.org/licenses/MIT):
//
// Copyright 2018-2019, CWI, TU Munich
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files
// (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// - The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
// IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// You can contact the authors via the FSST source repository : https://github.com/cwida/fsst
#include "perfevent/PerfEvent.hpp"
#include "fsst/fsst.h"
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <lz4.h>

using namespace std;

enum class AlgType : uint8_t {
  cpp_find,
  naive_kmp_on_compressed_data
};

std::string type_to_string(AlgType alg) {
  switch (alg) {
    case AlgType::cpp_find:
      return "c++-find";
    case AlgType::naive_kmp_on_compressed_data:
      return "naive-kmp[compressed]";
    default:
      return "Unknown";
  }
}

/// Base class for all compression tests.
class CompressionRunner {
   public:
   /// Store the compressed corpus. Returns the compressed size
   virtual uint64_t compressCorpus(const vector<string>& data, unsigned long &bareSize, double &bulkTime, double& compressionTime, bool verbose) = 0;
   /// Decompress some selected rows, separated by newlines. The line number are in ascending order. The target buffer is guaranteed to be large enough
   virtual uint64_t decompressRows(vector<char>& target, const vector<unsigned>& lines) = 0;
   /// Run Like.
   virtual uint64_t runLike(std::vector<char>& target, const std::string& pattern, AlgType algType) = 0;
};

/// No compresssion. Just used for debugging
class NoCompressionRunner : public CompressionRunner {
   private:
   /// The uncompressed data
   vector<string> data;

   public:
   /// Store the compressed corpus. Returns the compressed size
   uint64_t compressCorpus(const vector<string>& data, unsigned long& bareSize, double& bulkTime, double& compressionTime, bool /*verbose*/) override {
      auto startTime = std::chrono::steady_clock::now();
      this->data = data;
      uint64_t result = sizeof(uint32_t);
      for (auto& d : data)
         result += d.length() + sizeof(uint32_t);
      auto stopTime = std::chrono::steady_clock::now();
      bareSize = result;
      bulkTime = compressionTime = std::chrono::duration<double>(stopTime - startTime).count();
      return result;
   }
   /// Decompress some selected rows, separated by newlines. The line number are in ascending order. The target buffer is guaranteed to be large enough
   virtual uint64_t decompressRows(vector<char>& target, const vector<unsigned>& lines) {
      char* writer = target.data();
      for (auto l : lines) {
         auto& s = data[l];
         auto len = s.length();
         memcpy(writer, s.data(), len);
         writer[len] = '\n';
         writer += len + 1;
      }
      return writer - target.data();
   }

  uint64_t run_cpp_find(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    for (unsigned index = 0, limit = data.size(); index != limit; ++index) {
      if (data[index].find(pattern) != std::string::npos) {
        auto len = data[index].size();
        memcpy(writer, data[index].data(), len);
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return writer - target.data();
  }

  uint64_t runLike(std::vector<char>& target, const std::string& pattern, AlgType algType) override {
    switch (algType) {
      case AlgType::cpp_find: return run_cpp_find(target, pattern);
      default: return 0;
    }
  }
};

/// FSST compression
class FSSTCompressionRunner : public CompressionRunner {
   private:
   /// The decode
   fsst_decoder_t decoder;
   /// The compressed data
   vector<unsigned char> compressedData;
   /// The offsets
   vector<unsigned> offsets;

   public:
   FSSTCompressionRunner() {}
   FSSTCompressionRunner(unsigned /*blockSizeIgnored*/) {}

   /// Store the compressed corpus. Returns the compressed size
   uint64_t compressCorpus(const vector<string>& data, unsigned long& bareSize, double& bulkTime, double& compressionTime, bool verbose) override {
      compressedData.clear();
      offsets.clear();

      vector<unsigned long> rowLens, compressedRowLens;
      vector<unsigned char*> rowPtrs, compressedRowPtrs;
      rowLens.reserve(data.size());
      compressedRowLens.resize(data.size());
      rowPtrs.reserve(data.size());
      compressedRowPtrs.resize(data.size() + 1);
      unsigned long totalLen = 0;
      for (auto& d : data) {
         totalLen += d.size();
         rowLens.push_back(d.size());
         rowPtrs.push_back(reinterpret_cast<unsigned char*>(const_cast<char*>(d.data())));
      }

      auto firstTime = std::chrono::steady_clock::now();
      vector<unsigned long> dummy;
      if (getenv("LOOP"))
         for (int i = 0; i < 10000; i++) fsst_destroy(fsst_create(data.size(), rowLens.data(), const_cast<const unsigned char**>(rowPtrs.data()), false));
      auto encoder = fsst_create(data.size(), rowLens.data(), const_cast<const unsigned char**>(rowPtrs.data()), false);
      auto createTime = std::chrono::steady_clock::now();
      vector<unsigned char> compressionBuffer, fullBuffer;
      fullBuffer.resize(totalLen);
      unsigned char *fullBuf = fullBuffer.data();
      unsigned stringEnd = 0;
      for (auto& d : data) {
         memcpy(fullBuf + stringEnd, d.data(), d.length());
         stringEnd += d.length();
      }
      compressionBuffer.resize(16 + 2 * totalLen);
      auto copyTime = std::chrono::steady_clock::now();
      const unsigned char* fullBufPtr = fullBuf;
      fsst_compress(encoder, 1, &totalLen, &fullBufPtr, compressionBuffer.size(), compressionBuffer.data(), compressedRowLens.data(), compressedRowPtrs.data());
      auto startTime = std::chrono::steady_clock::now();
      fsst_compress(encoder, data.size(), rowLens.data(), const_cast<const unsigned char**>(rowPtrs.data()), compressionBuffer.size(), compressionBuffer.data(), compressedRowLens.data(), compressedRowPtrs.data());
      auto stopTime = std::chrono::steady_clock::now();
      unsigned long compressedLen = data.empty() ? 0 : (compressedRowPtrs[data.size() - 1] + compressedRowLens[data.size() - 1] - compressionBuffer.data());

      compressedData.resize(compressedLen + 8192);
      memcpy(compressedData.data(), compressionBuffer.data(), compressedLen);
      offsets.reserve(data.size());
      compressedRowPtrs[data.size()] = compressionBuffer.data() + compressedLen;
      for (unsigned index = 0, limit = data.size(); index != limit; ++index)
         offsets.push_back(compressedRowPtrs[index + 1] - compressionBuffer.data());
      bareSize = compressedData.size();
      uint64_t result = bareSize + (offsets.size() * sizeof(unsigned));
      {
         unsigned char buffer[sizeof(fsst_decoder_t)];
         unsigned dictLen = fsst_export(encoder, buffer);
         fsst_destroy(encoder);
         result += dictLen;

         fsst_import(&decoder, buffer);
      }
      double oneTime = std::chrono::duration<double>(createTime - firstTime).count();
      bulkTime = std::chrono::duration<double>(startTime - copyTime).count();
      compressionTime = std::chrono::duration<double>(stopTime - startTime).count();
      if (verbose) {
         cout << "# symbol table construction time: " << oneTime << endl;
         cout << "# compress-bulk time: " << bulkTime << endl;
         cout << "# compress time: " << compressionTime << endl;
      }
      bulkTime += oneTime;
      compressionTime += oneTime;

      return result;
   }
   /// Decompress some selected rows, separated by newlines. The line number are in ascending order. The target buffer is guaranteed to be large enough
   virtual uint64_t decompressRows(vector<char>& target, const vector<unsigned>& lines) {
      char* writer = target.data();
      auto limit = writer + target.size();

      // std::cerr << "target.size=" << target.size() << std::endl;

      auto data = compressedData.data();
      auto offsets = this->offsets.data();
      for (auto l : lines) {
         auto start = l ? offsets[l - 1] : 0, end = offsets[l];
         unsigned len = fsst_decompress(&decoder, end - start, data + start, limit - writer, reinterpret_cast<unsigned char*>(writer));
         writer[len] = '\n';
         writer += len + 1;
      }
      return writer - target.data();
   }

  uint64_t runLike(vector<char>& target, const std::string& pattern, AlgType alg_type) override {
    return 0;
  }
};

static pair<bool, tuple<double, double>> doFullDecompression(CompressionRunner& runner, string file, const unsigned num_repeats=100, bool verbose=false) {
  uint64_t totalSize = 0;
  bool debug = getenv("DEBUG");
  NoCompressionRunner debugRunner;

  // Read the corpus
  vector<string> corpus;
  uint64_t corpusLen = 0;
  {
    ifstream in(file);
    if (!in.is_open()) {
        cerr << "unable to open " << file << endl;
        return {false, {}};
    }
    string line;
    while (getline(in, line)) {
        corpusLen += line.length() + 1;
        corpus.push_back(move(line));
        if (corpusLen > 7000000) break;
    }
  }
  corpusLen += 4096;

  // Compress it
  double bulkTime, compressionTime;
  unsigned long bareSize;
  totalSize += runner.compressCorpus(corpus, bareSize, bulkTime, compressionTime, verbose);
  if (debug) {
    double ignored;
    debugRunner.compressCorpus(corpus, bareSize, ignored, ignored, false);
  }

  // Test different selectivities
  vector<char> targetBuffer, debugBuffer;
  targetBuffer.resize(corpusLen);
  if (debug) debugBuffer.resize(corpusLen);

  vector<unsigned> row_indices;
  for (unsigned index = 0, limit = corpus.size(); index != limit; ++index)
      row_indices.push_back(index);

  for (unsigned index = 0; index != num_repeats; ++index)
    runner.decompressRows(targetBuffer, row_indices);

  auto startTime = std::chrono::steady_clock::now();
  for (unsigned index = 0; index != num_repeats; ++index)
    runner.decompressRows(targetBuffer, row_indices);
  auto stopTime = std::chrono::steady_clock::now();

  auto timing_info = tuple<double, unsigned>(
    std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime).count(),
    corpusLen
  );

  if (verbose)
    std::cout << "# total compress size: " << totalSize << endl;
  tuple<double, double> result;

  auto time = std::get<0>(timing_info);
  auto size = std::get<1>(timing_info);
  return {true, {
    time / num_repeats / 1'000'000,
    (size / time) * num_repeats * 1'000'000'000
  }};
}

static pair<bool, tuple<double, double>> doLike(CompressionRunner& runner, string file, const string pattern, AlgType algType, const unsigned num_repeats=100, bool verbose=false) {
  uint64_t totalSize = 0;
  bool debug = getenv("DEBUG");
  NoCompressionRunner debugRunner;

  // Read the corpus
  vector<string> corpus;
  uint64_t corpusLen = 0;
  {
    ifstream in(file);
    if (!in.is_open()) {
        cerr << "unable to open " << file << endl;
        return {false, {}};
    }
    string line;
    while (getline(in, line)) {
        corpusLen += line.length() + 1;
        corpus.push_back(move(line));
        if (corpusLen > 7000000) break;
    }
  }
  corpusLen += 4096;

  // Compress it
  double bulkTime, compressionTime;
  unsigned long bareSize;
  totalSize += runner.compressCorpus(corpus, bareSize, bulkTime, compressionTime, verbose);
  if (debug) {
    double ignored;
    debugRunner.compressCorpus(corpus, bareSize, ignored, ignored, false);
  }

  // Test different selectivities
  vector<char> targetBuffer, debugBuffer;
  targetBuffer.resize(corpusLen);
  if (debug) debugBuffer.resize(corpusLen);

  for (unsigned index = 0; index != num_repeats; ++index)
    runner.runLike(targetBuffer, pattern, algType);

  auto startTime = std::chrono::steady_clock::now();
  for (unsigned index = 0; index != num_repeats; ++index)
    runner.runLike(targetBuffer, pattern, algType);
  auto stopTime = std::chrono::steady_clock::now();

  auto timing_info = tuple<double, unsigned>(
    std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime).count(),
    corpusLen
  );

  if (verbose)
    std::cout << "# total compress size: " << totalSize << endl;
  tuple<double, double> result;

  auto time = std::get<0>(timing_info);
  auto size = std::get<1>(timing_info);
  return {true, {
    time / num_repeats / 1'000'000,
    (size / time) * num_repeats * 1'000'000'000
  }};
}

int main(int argc, const char* argv[]) {
   if (argc < 3)
      return -1;

   string method = argv[1];
   std::string pattern = argv[2];
   vector<string> files;
   for (int index = 3; index < argc; ++index) {
      string f = argv[index];
      if (f == "--exclude") {
         auto iter = find(files.begin(), files.end(), argv[++index]);
         if (iter != files.end()) files.erase(iter);
      } else {
         files.push_back(move(f));
      }
   }

  if (method == "full-decompression") {
    assert(files.size() == 1);

    NoCompressionRunner runner1;
    auto r1 = doFullDecompression(runner1, files.front());
    assert(r1.first);

    FSSTCompressionRunner runner2;
    auto r2 = doFullDecompression(runner2, files.front());
    assert(r2.first);

    cout << "type \t time [ms] \t throughput [#tuples / s]" << endl;
    cout << "vanilla" << "\t" << std::get<0>(r1.second) << "\t" << std::get<1>(r1.second) << std::endl;
    cout << "fsst" << "\t" << std::get<0>(r2.second) << "\t" << std::get<1>(r2.second) << std::endl;
  } else if (method == "like") {
    // NoCompressionRunner runner;
    // assert(files.size() == 1);
    // std::cerr << "file=" << files.front() << std::endl;
    // auto r1 = doLike(runner, files.front(), pattern, AlgType::cpp_find); 
    // assert(r1.first);

    // NoCompressionRunner runner;
    // assert(files.size() == 1);
    // std::cerr << "file=" << files.front() << std::endl;
    // auto r1 = doLike(runner, files.front(), pattern, AlgType::cpp_find); 
    // assert(r1.first);

    // cout << "algo\tthroughput [#tuples / s]" << endl;
    // cout << type_to_string(AlgType::cpp_find) << "\t" << std::get<0>(r.second) << "\t" << std::get<1>(r.second) << std::endl;
  } else {
    cerr << "unknown method " << method << endl;
    return 1;
  }
}
