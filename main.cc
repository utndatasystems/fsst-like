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
#include <regex>
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
  cpp_regex,
  cpp_memmem,
  kmp_on_decompressed_data,
  kmp_on_compressed_data
};

std::string type_to_string(AlgType alg) {
  switch (alg) {
    case AlgType::cpp_find:
      return "c++-find";
    case AlgType::cpp_regex:
      return "c++-regex";
    case AlgType::cpp_memmem:
      return "c++-memmem";
    case AlgType::kmp_on_decompressed_data:
      return "kmp[decompressed]";
    case AlgType::kmp_on_compressed_data:
      return "kmp[compressed]";
    default:
      return "Unknown";
  }
}

static constexpr unsigned FSST_SIZE = 255;
#define FSST_CORRUPT 32774747032022883 /* 7-byte number in little endian containing "corrupt" */

class StateMachine {
public:
  StateMachine(const std::string& pattern) : P(pattern), m(pattern.size()) {
    build_pi();
  }

  void init_state_lookup(fsst_decoder_t& fsst_decoder) {
    // Resize.
    fsst_symbols.reserve(FSST_SIZE);

    // Take the FSST symbols.
    for (unsigned index = 0; index < FSST_SIZE; ++index) {
      auto len = fsst_decoder.len[index];
      auto raw_symbol = fsst_decoder.symbol[index];

      if (raw_symbol == FSST_CORRUPT) {
        continue;
      }

      std::string symbol;
      for (unsigned k = 0; k != len; ++k) {
	  		int c = (raw_symbol >> 8*k) & 0xFF;
        symbol += static_cast<char>(c);
      }

      // Add symbol.
      fsst_symbols.push_back(symbol);
    }
    
    // Init the fsst_size.
    fsst_size = fsst_symbols.size();

    // Init the state lookup.
    state_lookup.resize(P.size() * fss);
  }

  void build_state_lookup(fsst_decoder_t& fsst_decoder) {
    // Init the state cache.
    init_state_lookup(fsst_decoder);

    for (unsigned index = 0; index != m; ++index) {
      for (unsigned symbol_index = 0; symbol_index != fsst_size; ++index) {
        init_state(index);

        accept_symbol(symbol);

        state_lookup[i * fsst_size + symbol_index] = curr_state;
      }
    }
  }

    for i in range(len(self.pattern)):
      for symbol in fsst_symbols:
        # Set the current state to position `i`.
        self.init_state(i)

        # Simulate.
        self.accept_symbol(symbol)

        # Cache the state we arrived at.
        self.cache_state[i * len(fsst_symbols) + symbol.index] = self.curr_state
        
  void init_state(unsigned pos = 0) {
    curr_state = pos;
  }

  void accept(char c) {
    while ((curr_state > 0) && (P[curr_state] != c)) {
      curr_state = pi[curr_state - 1];
    }

    if (P[curr_state] == c)
      ++curr_state;
  }

  void accept(unsigned symbol_index) {
    while ((curr_state > 0) && (P[curr_state] != c)) {
      curr_state = pi[curr_state - 1];
    }

    if (P[curr_state] == c)
      ++curr_state;
  }

  // void accept_symbol(symbol) {
  //   for letter in symbol.val:
  //     self.accept_letter(letter)

  //     # Already reached the final state?
  //     if self.curr_state == len(self.pattern):
  //       return
      

  bool match(const std::string& text) {
    // Init.
    init_state();

    // Iterate.
    for (unsigned index = 0, limit = text.size(); index != limit; ++index) {
      accept(text[index]);

      if (curr_state == m)
        return true;
    }
    return false;
  }

  bool match(const char* ptr, unsigned len) {
    // Init.
    init_state();

    // Iterate.
    for (unsigned index = 0, limit = len; index != limit; ++index) {
      accept(ptr[index]);

      if (curr_state == m)
        return true;
    }
    return false;
  }

private:
  std::string P;
  unsigned m;
  unsigned curr_state;
  std::vector<unsigned> pi;
  unsigned fsst_size;
  std::vector<std::string> fsst_symbols;
  std::vector<unsigned> state_lookup;

  void build_pi() {
    pi.assign(P.size(), 0);

    pi[0] = 0;

    unsigned k = 0;
    for (unsigned q = 1, limit = P.size(); q != limit; ++q) {
      while ((k > 0) && (P[k] != P[q])) {
        k = pi[k - 1];
      }

      // Match? Then advance.
      if (P[k] == P[q])
        ++k;
      
      // Store the state.
      pi[q] = k;
    }
  }
};

/// Base class for all compression tests.
class CompressionRunner {
   public:
   /// Store the compressed corpus. Returns the compressed size
   virtual uint64_t compressCorpus(const vector<string>& data, unsigned long &bareSize, double &bulkTime, double& compressionTime, bool verbose) = 0;
   /// Decompress some selected rows, separated by newlines. The line number are in ascending order. The target buffer is guaranteed to be large enough
   virtual uint64_t decompressRows(vector<char>& target, const vector<unsigned>& lines) = 0;
   /// Run `LIKE`.
   virtual size_t runLike(std::vector<char>& target, const std::string& pattern, AlgType algType) = 0;
};

/// No compresssion. Just used for debugging
class NoCompressionRunner : public CompressionRunner {
   private:
   /// The uncompressed data
   vector<string> data;

   public:
   /// Store the compressed corpus. Returns the compressed size
   uint64_t compressCorpus(const vector<string>& data, unsigned long& bareSize, double& bulkTime, double& compressionTime, bool /*verbose*/) override {
      auto startTime = std::chrono::high_resolution_clock::now();
      this->data = data;
      uint64_t result = sizeof(uint32_t);
      for (auto& d : data)
         result += d.length() + sizeof(uint32_t);
      auto stopTime = std::chrono::high_resolution_clock::now();
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

  size_t run_cpp_find(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    size_t count = 0;
    for (unsigned index = 0, limit = data.size(); index != limit; ++index) {
      if (data[index].find(pattern) != std::string::npos) {
        ++count;
        auto len = data[index].size();
        memcpy(writer, data[index].data(), len);
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  uint64_t run_cpp_regex(std::vector<char>& target, const std::string& pattern) {
    std::regex regex_pattern(pattern, std::regex_constants::optimize);
    char* writer = target.data();

    size_t count = 0;
    for (unsigned index = 0, limit = data.size(); index != limit; ++index) {
      if (std::regex_search(data[index], regex_pattern)) {
        auto len = data[index].size();
        memcpy(writer, data[index].data(), len);
        writer[len] = '\n';
        writer += len + 1;
      }
    }

    return count;
  }

  size_t run_kmp_on_decompressed_data(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();

    auto stateMachine = StateMachine(pattern);
    size_t count = 0;
    for (unsigned index = 0, limit = data.size(); index != limit; ++index) {
       if (stateMachine.match(data[index])) {
        ++count;
        auto len = data[index].size();
        memcpy(writer, data[index].data(), len);
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  size_t run_cpp_memmem(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    size_t count = 0;
    size_t pattern_len = pattern.size();
    for (unsigned index = 0, limit = data.size(); index != limit; ++index) {
      if (memmem(data[index].data(), data[index].size(), pattern.data(), pattern_len)) {
        ++count;
        auto len = data[index].size();
        memcpy(writer, data[index].data(), len);
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  uint64_t runLike(std::vector<char>& target, const std::string& pattern, AlgType algType) override {
    switch (algType) {
      case AlgType::cpp_find: return run_cpp_find(target, pattern);
      case AlgType::cpp_regex: return run_cpp_regex(target, pattern);
      case AlgType::cpp_memmem: return run_cpp_memmem(target, pattern);
      case AlgType::kmp_on_decompressed_data: return run_kmp_on_decompressed_data(target, pattern);
      case AlgType::kmp_on_compressed_data: return 0;
      default: return 0;
    }
  }
};

/// FSST compression
class FSSTCompressionRunner : public CompressionRunner {
   private:
   /// The data size
   std::size_t data_size;
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

      // Set the data size.
      data_size = data.size();

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

      auto firstTime = std::chrono::high_resolution_clock::now();
      vector<unsigned long> dummy;
      if (getenv("LOOP"))
         for (int i = 0; i < 10000; i++) fsst_destroy(fsst_create(data.size(), rowLens.data(), const_cast<const unsigned char**>(rowPtrs.data()), false));
      auto encoder = fsst_create(data.size(), rowLens.data(), const_cast<const unsigned char**>(rowPtrs.data()), false);
      auto createTime = std::chrono::high_resolution_clock::now();
      vector<unsigned char> compressionBuffer, fullBuffer;
      fullBuffer.resize(totalLen);
      unsigned char *fullBuf = fullBuffer.data();
      unsigned stringEnd = 0;
      for (auto& d : data) {
         memcpy(fullBuf + stringEnd, d.data(), d.length());
         stringEnd += d.length();
      }
      compressionBuffer.resize(16 + 2 * totalLen);
      auto copyTime = std::chrono::high_resolution_clock::now();
      const unsigned char* fullBufPtr = fullBuf;
      fsst_compress(encoder, 1, &totalLen, &fullBufPtr, compressionBuffer.size(), compressionBuffer.data(), compressedRowLens.data(), compressedRowPtrs.data());
      auto startTime = std::chrono::high_resolution_clock::now();
      fsst_compress(encoder, data.size(), rowLens.data(), const_cast<const unsigned char**>(rowPtrs.data()), compressionBuffer.size(), compressionBuffer.data(), compressedRowLens.data(), compressedRowPtrs.data());
      auto stopTime = std::chrono::high_resolution_clock::now();
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
      auto writer_limit = writer + target.size();

      // std::cerr << "target.size=" << target.size() << std::endl;

      auto data = compressedData.data();
      auto offsets = this->offsets.data();
      for (auto l : lines) {
         auto start = l ? offsets[l - 1] : 0, end = offsets[l];
         unsigned len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));
         writer[len] = '\n';
         writer += len + 1;
      }
      return writer - target.data();
   }

  size_t run_cpp_find(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];
      unsigned len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));
      
      // Match?
      if (std::string_view(writer, len).find(pattern) != std::string::npos) {
        ++count;
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  size_t run_cpp_regex(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    std::regex regex_pattern(pattern, std::regex_constants::optimize);

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];
      unsigned len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));

      // Match?
      if (std::regex_search(writer, writer + len, regex_pattern)) {
        ++count;
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  size_t run_cpp_memmem(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    size_t pattern_len = pattern.size();
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];
      unsigned len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));
      
      // Match?
      if (memmem(writer, len, pattern.data(), pattern_len)) {
        ++count;
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  size_t run_kmp_on_decompressed_data(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    auto stateMachine = StateMachine(pattern);

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];
      unsigned len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));

      // Match?
      if (stateMachine.match(writer, len)) {
        ++count;
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  size_t run_kmp_on_compressed_data(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    auto stateMachine = StateMachine(pattern);

    stateMachine.init_state_lookup(decoder);

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];
      unsigned len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));

      // Match?
      if (stateMachine.match(writer, len)) {
        ++count;
        writer[len] = '\n';
        writer += len + 1;
      }
    }
    return count;
  }

  uint64_t runLike(vector<char>& target, const std::string& pattern, AlgType alg_type) override {
    switch (alg_type) {
      case AlgType::cpp_find: return run_cpp_find(target, pattern);
      case AlgType::cpp_regex: return run_cpp_regex(target, pattern);
      case AlgType::cpp_memmem: return run_cpp_memmem(target, pattern);
      case AlgType::kmp_on_decompressed_data: return run_kmp_on_decompressed_data(target, pattern);
      case AlgType::kmp_on_compressed_data: return run_kmp_on_compressed_data(target, pattern);
      default: return 0;
    }
  }
};

static pair<bool, tuple<size_t, double, size_t, unsigned>> doFullDecompression(CompressionRunner& runner, string file, const unsigned num_repeats=100, bool verbose=false) {
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
        // if (corpusLen > 7000000) break;
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

  auto startTime = std::chrono::high_resolution_clock::now();
  for (unsigned index = 0; index != num_repeats; ++index)
    runner.decompressRows(targetBuffer, row_indices);
  auto stopTime = std::chrono::high_resolution_clock::now();

  return {true, {
    0,
    std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime).count(),
    corpus.size(),
    num_repeats
  }};
}

static pair<bool, tuple<size_t, double, size_t, unsigned>> doLike(CompressionRunner& runner, string file, const string pattern, AlgType algType, const unsigned num_repeats=1, bool verbose=false) {
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

  auto count = 0;
  for (unsigned index = 0; index != num_repeats; ++index)
    count = runner.runLike(targetBuffer, pattern, algType);

  auto startTime = std::chrono::high_resolution_clock::now();
  for (unsigned index = 0; index != num_repeats; ++index)
    count = runner.runLike(targetBuffer, pattern, algType);
  auto stopTime = std::chrono::high_resolution_clock::now();

  return {true, {
    count,
    std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime).count(),
    corpus.size(),
    num_repeats
  }};
}

std::tuple<size_t, double, double, unsigned> display(std::pair<bool, std::tuple<size_t, double, size_t, unsigned>> r) {
  auto count = std::get<0>(r.second);
  auto time = std::get<1>(r.second);
  auto size = std::get<2>(r.second);
  auto num_repeats = std::get<3>(r.second);
  return {
    count,
    time / num_repeats / 1'000'000,
    1'000'000'000 * num_repeats * size / time,
    num_repeats
  };
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

    auto info1 = display(r1);
    auto info2 = display(r2);

    // Same count.
    assert(std::get<0>(info1) == std::get<0>(info2));

    cout << "type\ttime [ms]\t throughput [#tuples / s]" << endl;
    cout << "vanilla" << "\t" << std::get<1>(info1) << "\t" << std::get<2>(info1) << std::endl;
    cout << "fsst" << "\t" << std::get<1>(info2) << "\t" << std::get<2>(info2) << std::endl;
  } else if (method == "like") {
    assert(files.size() == 1);

    cout << "algo\ttype\ttime [ms]\t throughput [#tuples / s]" << endl;
    for (auto algType : {
      AlgType::cpp_find,
      AlgType::cpp_memmem,
      AlgType::kmp_on_decompressed_data,
      AlgType::kmp_on_compressed_data
    }) {
      NoCompressionRunner runner1;
      auto r1 = doLike(runner1, files.front(), pattern, algType);
      assert(r1.first);

      FSSTCompressionRunner runner2;
      auto r2 = doLike(runner2, files.front(), pattern, algType);
      assert(r2.first);

      auto info1 = display(r1);
      auto info2 = display(r2);

      // Same count.
      assert(std::get<0>(info1) == std::get<0>(info2));

      auto algo = type_to_string(algType);
      cout << algo << "\t" << "vanilla" << "\t" << std::get<1>(info1) << "\t" << std::get<2>(info1) << std::endl;
      cout << algo << "\t" << "fsst" << "\t" << std::get<1>(info2) << "\t" << std::get<2>(info2) << std::endl;
    }
  } else {
    cerr << "unknown method " << method << endl;
    return 1;
  }
}
