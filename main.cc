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
#include "Memory.hpp"
#include <cassert>
#include <regex>
#include <limits>
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
  kmp_lower_bound,
  kmp_on_decompressed_data,
  kmp_on_compressed_data,
  lookup_kmp_on_compressed_data,
  simple_simd_find,
};

std::string alg_type_to_string(AlgType alg) {
  switch (alg) {
    case AlgType::cpp_find:
      return "c++-find";
    case AlgType::cpp_regex:
      return "c++-regex";
    case AlgType::cpp_memmem:
      return "c++-memmem";
    case AlgType::kmp_lower_bound:
      return "lb-kmp[compressed]";
    case AlgType::kmp_on_decompressed_data:
      return "kmp[decompressed]";
    case AlgType::kmp_on_compressed_data:
      return "kmp[compressed]";
    case AlgType::lookup_kmp_on_compressed_data:
      return "lookup-kmp[compressed]";
    case AlgType::simple_simd_find:
      return "simple_simd_find";
    default:
      return "Unknown";
  }
}

static constexpr size_t infty = std::numeric_limits<size_t>::max();
static constexpr unsigned FSST_SIZE = 255;
#define FSST_CORRUPT 32774747032022883 /* 7-byte number in little endian containing "corrupt" */

/* Decompress a single string, inlined for speed. */
template<typename ConsumeCode, typename ConsumeChar>
inline bool  /* OUT: bytesize of the decompressed string. If > size, the decoded output is truncated to size. */
fsst_iterate(
  const fsst_decoder_t *decoder,  /* IN: use this symbol table for compression. */
  size_t lenIn,             /* IN: byte-length of compressed string. */
  const unsigned char *strIn,     /* IN: compressed string. */
  size_t size,              /* IN: byte-length of output buffer. */
  ConsumeChar&& consume_char,
  ConsumeCode&& consume_code
) {
   unsigned char*__restrict__ len = (unsigned char* __restrict__) decoder->len;
  //  unsigned char*__restrict__ strOut = (unsigned char* __restrict__) output;
   unsigned long long*__restrict__ symbol = (unsigned long long* __restrict__) decoder->symbol; 
   size_t code, posOut = 0, posIn = 0;
#ifndef FSST_MUST_ALIGN /* defining on platforms that require aligned memory access may help their performance */
#define FSST_UNALIGNED_STORE(dst,src) memcpy((unsigned long long*) (dst), &(src), sizeof(unsigned long long))
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)

  // TODO: Remove the stuff with `posOut`. Like everything.
   while (posOut+32 <= size && posIn+4 <= lenIn) {
      unsigned int nextBlock, escapeMask;
      memcpy(&nextBlock, strIn+posIn, sizeof(unsigned int));
      escapeMask = (nextBlock&0x80808080u)&((((~nextBlock)&0x7F7F7F7Fu)+0x7F7F7F7Fu)^0x80808080u);
      if (escapeMask == 0) {
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; //FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; //FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
     } else { 
         unsigned long firstEscapePos=__builtin_ctzl((unsigned long long) escapeMask)>>3;
         switch(firstEscapePos) { /* Duff's device */
         case 3: code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code];
                 // fall through
         case 2: code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code];
                 // fall through
         case 1: code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code];
                 // fall through
         case 0: posIn+=2; if (!consume_char(strIn[posIn - 1])) return true; posOut++; // strOut[posOut++] =  strIn[posIn-1]; /* decompress an escaped byte */
         }
      }
   }
   if (posOut+32 <= size) { // handle the possibly 3 last bytes without a loop
      if (posIn+2 <= lenIn) {
        if (strIn[posIn] != FSST_ESC) {
            code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            if (strIn[posIn] != FSST_ESC) {
               code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
            } else { 
               posIn += 2; if (!consume_char(strIn[posIn - 1])) return true; posOut++; // strOut[posOut++] = strIn[posIn-1]; 
            }
         } else {
            // Consume the regualr byte.
            if (!consume_char(strIn[posIn + 1]))
              return true;
            
            posIn += 2; posOut++; // posOut++;
         } 
      }
      if (posIn < lenIn) { // last code cannot be an escape
         code = strIn[posIn++]; if (!consume_code(code)) return true; posOut += len[code]; // FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); posOut += len[code]; 
      }
   }
#else
   while (posOut+8 <= size && posIn < lenIn)
      if ((code = strIn[posIn++]) < FSST_ESC) { /* symbol compressed as code? */
         FSST_UNALIGNED_STORE(strOut+posOut, symbol[code]); /* unaligned memory write */
         posOut += len[code];
      } else { 
         strOut[posOut] = strIn[posIn]; /* decompress an escaped byte */
         posIn++; posOut++; 
      }
#endif
#endif
   while (posIn < lenIn)
      if ((code = strIn[posIn++]) < FSST_ESC) {
         size_t posWrite = posOut, endWrite = posOut + len[code];
         unsigned char* __restrict__ symbolPointer = ((unsigned char* __restrict__) &symbol[code]) - posWrite;
         if ((posOut = endWrite) > size) endWrite = size;
         for(; posWrite < endWrite; posWrite++) { /* only write if there is room */
            if (!consume_char(symbolPointer[posWrite])) return true; // strOut[posWrite] = symbolPointer[posWrite];
         }
      } else {
         if (posOut < size) {
            if (!consume_char(strIn[posIn])) return true; // strOut[posOut] = strIn[posIn]; /* idem */
         }
         posIn++; posOut++; 
      } 
   if (posOut >= size && (decoder->zeroTerminated&1)) {
     assert(0); //  strOut[size-1] = 0;
   }
  //  std::cerr << std::endl;
   return false; /* full size of decompressed string (could be >size, then the actually decompressed part) */
}

class StateMachine {
public:
  StateMachine(const std::string& pattern) : P(pattern), m(pattern.size()) {
    build_pi();
  }

  void init_fsst_symbols(fsst_decoder_t& fsst_decoder) {
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
  }

  void build_lookup_table() {
    // Init the lookup table.
    lookup_table.resize(P.size() * fsst_size);

    for (unsigned index = 0; index != m; ++index) {
      for (unsigned code = 0; code != fsst_size; ++code) {
        // Init the state.
        init_state(index);

        accept_symbol(code);

        lookup_table[index * fsst_size + code] = curr_state;
      }
    }
  }

  void init_state(unsigned pos = 0) {
    curr_state = pos;
  }

  void accept(char c) {
    // TODO: Maybe optimize this when `curr_state` is anyway always 0.
    // TODO: Like we should optimize for the last if.
    while ((curr_state > 0) && (P[curr_state] != c)) {
      curr_state = pi[curr_state - 1];
    }

    if (P[curr_state] == c)
      ++curr_state;
  }

  void accept_symbol(size_t code) {
    for (auto c : fsst_symbols[code]) {
      accept(c);

      if (curr_state == m)
        return;
    }
  }

  void accept_symbol_with_lookup(size_t code) {
    // Use the state lookup table.
    curr_state = lookup_table[curr_state * fsst_size + code];
  }

  bool match(const char* ptr, unsigned len) {
    // Init.
    init_state();

    // Iterate.
    for (unsigned index = 0, limit = len; index != limit; ++index) {
      accept(ptr[index]);

      // TODO: Here, we can have two paths.
      // TODO: If `curr_state` is very far from the end, don't have this `if` here.
      // TODO: Otherwise, have it.
      // TODO: But maybe only when necessary (like traverse until some point and then switch).
      if (curr_state == m)
        return true;
    }
    return false;
  }

  // Another implementation.
  bool match(const std::string& text) { return match(text.data(), text.size()); }

  bool inline_match(const fsst_decoder_t* decoder, size_t lenIn, const unsigned char* strIn, size_t size) {
    // Init.
    init_state();

    auto consume_char = [this](unsigned char c) {
      accept(c);
      return curr_state != m;
    };

    auto consume_code = [this](size_t code) {
      // TODO: Remove this assert at the end.
      assert(code < fsst_size);
      accept_symbol(code);
      return curr_state != m;
    };

    return fsst_iterate(decoder, lenIn, strIn, size, consume_char, consume_code);
  }

  bool inline_match_with_lookup(const fsst_decoder_t* decoder, size_t lenIn, const unsigned char* strIn, size_t size) {
    // Init.
    init_state();

    auto consume_char = [this](unsigned char c) {
      accept(c);
      return curr_state != m;
    };

    auto consume_code = [this](size_t code) {
      // TODO: Remove this assert at the end.
      assert(code < fsst_size);
      accept_symbol_with_lookup(code);
      return curr_state != m;
    };

    return fsst_iterate(decoder, lenIn, strIn, size, consume_char, consume_code);
  }

private:
  std::string P;
  unsigned m;
  unsigned curr_state;
  std::vector<unsigned> pi;
  unsigned fsst_size;
  std::vector<std::string> fsst_symbols;
  std::vector<unsigned> lookup_table;

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

class DummyStateMachine : public StateMachine {
public:
  using StateMachine::StateMachine;

  // This is the best we could achieve with KMP.
  bool inline_match(const fsst_decoder_t* decoder, size_t lenIn, const unsigned char* strIn, size_t size) {
    auto dummy_consume_char = [this](unsigned char c) {
      ++c;
      return false;
    };

    auto dummy_consume_code = [this](size_t code) {
      ++code;
      return false;
    };

    fsst_iterate(decoder, lenIn, strIn, size, dummy_consume_char, dummy_consume_code);
    return true;
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
   virtual size_t runLike(std::vector<char>& target, const std::string& pattern, AlgType algType, const std::vector<unsigned>& oracle) = 0;
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

  size_t run_simple_simd_find(std::vector<char>& target, const std::string& pattern) {
    if (pattern.size() > 8) {
      throw std::runtime_error("Pattern too long");
    }

    char* writer = target.data();
    size_t count = 0;
    size_t pattern_len = pattern.size();
    uint64_t pattern_block = Memory::LoadSafe8(pattern.data(), pattern.size());
    uint64_t pattern_mask = uint64_t(0xffffffffffffffff) >> (8 * (8 - pattern_len));
    for (unsigned index = 0, limit = data.size(); index != limit; ++index) {

      // Fast path for full 8 byte blocks in corps.
      const string& corps = data[index];
      unsigned corps_len = corps.size();
      unsigned i = 0;
      if(corps_len > 8) {
        for (; i < corps_len - 8; i++) {
          uint64_t corps_block = Memory::Load<uint64_t>(corps.data() + i);
          if (((corps_block ^ pattern_block) & pattern_mask) == 0) {
            goto match;
          }
        }
      }

      // Handle last 8 byte, individually.
      for (; i < corps_len - pattern_len + 1; i++) {
        uint64_t corps_block = Memory::LoadSafe8(corps.data() + i, corps_len - i);
        if (((corps_block ^ pattern_block) & pattern_mask) == 0) {
          goto match;
        }
      }

      continue;

      match:
      ++count;
      auto len = data[index].size();
      memcpy(writer, data[index].data(), len);
      writer[len] = '\n';
      writer += len + 1;
    }

    return count;
  }

  uint64_t runLike(std::vector<char>& target, const std::string& pattern, AlgType algType, const std::vector<unsigned>& /* oracle */) override {
    switch (algType) {
      case AlgType::cpp_find: return run_cpp_find(target, pattern);
      case AlgType::cpp_regex: return run_cpp_regex(target, pattern);
      case AlgType::cpp_memmem: return run_cpp_memmem(target, pattern);
      case AlgType::kmp_on_decompressed_data: return run_kmp_on_decompressed_data(target, pattern);
      case AlgType::simple_simd_find: return run_simple_simd_find(target, pattern);
      default: return infty;
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

    // Init the machine.
    auto stateMachine = StateMachine(pattern);
    
    // Init the FSST-symbols.
    stateMachine.init_fsst_symbols(decoder);

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];

      // Check.
      auto ans = stateMachine.inline_match(&decoder, end - start, data + start, writer_limit - writer);

      // Match?
      if (ans) {
        ++count;

        // Decompress.
        unsigned new_len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));
        writer[new_len] = '\n';
        writer += new_len + 1;
      }
    }
    return count;
  }

  size_t run_lookup_kmp_on_compressed_data(std::vector<char>& target, const std::string& pattern) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    // Init the machine.
    auto stateMachine = StateMachine(pattern);

    // Init the FSST symbols.
    stateMachine.init_fsst_symbols(decoder);

    // Build the state lookup table.
    stateMachine.build_lookup_table();

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];

      // Check.
      auto ans = stateMachine.inline_match_with_lookup(&decoder, end - start, data + start, writer_limit - writer);

      // Match?
      if (ans) {
        ++count;

        // Decompress.
        unsigned new_len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));
        writer[new_len] = '\n';
        writer += new_len + 1;
      }
    }
    return count;
  }

  size_t run_kmp_lower_bound(std::vector<char>& target, const std::string& pattern, const std::vector<unsigned>& oracle) {
    char* writer = target.data();
    auto writer_limit = writer + target.size();

    DummyStateMachine dummyStateMachine(pattern);

    auto data = compressedData.data();
    auto offsets = this->offsets.data();
    size_t count = 0;
    for (unsigned index = 0, limit = this->data_size; index != limit; ++index) {
      auto start = index ? offsets[index - 1] : 0, end = offsets[index];

      dummyStateMachine.inline_match(&decoder, end - start, data + start, writer_limit - writer);

      if ((count != oracle.size()) && (index == oracle[count])) {
        ++count;

        // Decompress.
        unsigned new_len = fsst_decompress(&decoder, end - start, data + start, writer_limit - writer, reinterpret_cast<unsigned char*>(writer));
        writer[new_len] = '\n';
        writer += new_len + 1;
      }
    }
    return count;
  }

  uint64_t runLike(vector<char>& target, const std::string& pattern, AlgType alg_type, const std::vector<unsigned>& oracle) override {
    switch (alg_type) {
      case AlgType::cpp_find: return run_cpp_find(target, pattern);
      case AlgType::cpp_regex: return run_cpp_regex(target, pattern);
      case AlgType::cpp_memmem: return run_cpp_memmem(target, pattern);
      case AlgType::kmp_lower_bound: return run_kmp_lower_bound(target, pattern, oracle);
      case AlgType::kmp_on_decompressed_data: return run_kmp_on_decompressed_data(target, pattern);
      case AlgType::kmp_on_compressed_data: return run_kmp_on_compressed_data(target, pattern);
      case AlgType::lookup_kmp_on_compressed_data: return run_lookup_kmp_on_compressed_data(target, pattern);
      default: return infty;
    }
  }
};

std::vector<unsigned> computeOracle(string file, const string pattern) {
  // Read the corpus
  std::vector<unsigned> oracle;
  {
    ifstream in(file);
    if (!in.is_open()) {
        cerr << "unable to open " << file << endl;
        return {false, {}};
    }
    string line;
    unsigned row_index = 0;
    while (getline(in, line)) {
      if (line.find(pattern) != std::string::npos) {
        oracle.push_back(row_index);
      }
      ++row_index;
    }
  }
  return oracle;
}

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

static pair<bool, tuple<size_t, double, size_t, unsigned>> doLike(CompressionRunner& runner, string file, const string pattern, AlgType algType, const std::vector<unsigned>& oracle, const unsigned num_repeats=1, bool verbose=false) {
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

  auto startTime = std::chrono::high_resolution_clock::now();
  auto count = 0;
  for (unsigned index = 0; index != num_repeats; ++index)
    count = runner.runLike(targetBuffer, pattern, algType, oracle);
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
    auto info_fsst = display(r2);

    // Same count.
    assert(std::get<0>(info1) == std::get<0>(info_fsst));

    cout << "type\ttime [ms]\t throughput [#tuples / s]" << endl;
    cout << "uncompressed" << "\t" << std::get<1>(info1) << "\t" << std::get<2>(info1) << std::endl;
    cout << "fsst" << "\t" << std::get<1>(info_fsst) << "\t" << std::get<2>(info_fsst) << std::endl;
  } else if (method == "like") {
    assert(files.size() == 1);

    // The ranking.
    std::vector<std::tuple<std::string, AlgType, double, double>> ranking;

    for (auto algType : {
      AlgType::simple_simd_find,
      AlgType::cpp_find,
      AlgType::cpp_memmem,
      AlgType::kmp_lower_bound,
      AlgType::kmp_on_decompressed_data,
      AlgType::kmp_on_compressed_data,
      AlgType::lookup_kmp_on_compressed_data,
    }) {
      std::cerr << "Running " << alg_type_to_string(algType) << ".." << std::endl;
      auto oracle = computeOracle(files.front(), pattern);

      NoCompressionRunner runner1;
      auto r1 = doLike(runner1, files.front(), pattern, algType, oracle);
      assert(r1.first);

      FSSTCompressionRunner runner2;
      auto r2 = doLike(runner2, files.front(), pattern, algType, oracle);
      assert(r2.first);

      auto info_uncompressed = display(r1);
      auto info_fsst = display(r2);

      std::cerr << std::get<0>(info_uncompressed) << " " << std::get<0>(info_fsst) << std::endl;

      // Check with the oracle.
      if (std::get<0>(info_uncompressed) != infty) assert(std::get<0>(info_uncompressed) == oracle.size());
      if (std::get<0>(info_fsst) != infty) assert(std::get<0>(info_fsst) == oracle.size());

      // Add uncompressed.
      if (std::get<0>(info_uncompressed) != infty) {
        ranking.push_back({
          "uncompressed",
          algType,
          std::get<1>(info_uncompressed),
          std::get<2>(info_uncompressed)
        });
      }

      // Add FSST.
      if (std::get<0>(info_fsst) != infty) {
        ranking.push_back({
          "fsst",
          algType,
          std::get<1>(info_fsst),
          std::get<2>(info_fsst)
        });
      }
    }

    // Sort by the time.
    std::sort(ranking.begin(), ranking.end(), [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
    });

    // And print.
    cout << "\nalgo, type, time [ms], throughput [#tuples / s]" << endl;
    for (auto elem : ranking) {
      auto data_type = std::get<0>(elem);
      auto algo = alg_type_to_string(std::get<1>(elem));
      cout << algo << ", " << data_type << ", " << std::get<2>(elem) << ", " << std::get<3>(elem) << std::endl;
    }
  } else {
    cerr << "unknown method " << method << endl;
    return 1;
  }
}
