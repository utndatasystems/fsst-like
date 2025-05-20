#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <set>
#include "BenchmarkDriver.hpp"
#include "StateMachine.hpp"
#include "Utility.hpp"
#include "src/SimdEverywhere.hpp"
// -------------------------------------------------------------------------------------
class SkippingEngine : public Engine {
public:
   SkippingEngine(std::string_view pattern, StateMachine2 state_machine)
       : pattern(pattern)
       , decode_buffer(128)
       , state_machine(std::move(state_machine))
   {
   }

   uint64_t skipped = 0;

   ~SkippingEngine()
   {
      std::cout << "skipped: " << skipped << std::endl;
      std::cout << "prepare time: " << (total / 1e6) << "ms" << std::endl;
   }

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result) final
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (block.GetRow(row_idx).find(pattern) != std::string_view::npos) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint64_t total = 0;

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result) final
   {
      state_machine.init_fsst_symbols(block.decoder);
      state_machine.build_lookup_table();

      auto begin = std::chrono::high_resolution_clock::now();
      std::vector<uint8_t> required_symbols = CreateRequiredSymbols(block);
      auto end = std::chrono::high_resolution_clock::now();
      total += std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

      if (required_symbols.size() <= 2) {
         // return SkippingScanPerRowSimple(block, result, required_symbols);
         return SkippingScanPerRowSimple(block, result, required_symbols);
      }
      // if (required_symbols.size() == 1) {
      //    return SkippingScanPerRowSimd1(block, result, required_symbols);
      // }
      // else if (required_symbols.size() == 2) {
      //    return SkippingScanPerRowSimd2(block, result, required_symbols);
      // }
      else {
         return NormalScan(block, result, required_symbols);
      }
   }

   uint32_t NormalScan(const FsstBlock& block, std::vector<uint32_t>& result, const std::vector<uint8_t>& required_symbols)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);

         // Matching code.
         const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(compressed_text.data());
         // bool match = state_machine.fsst_lookup_kmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         bool match = state_machine.fsst_lookup_zerokmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         if (match) {
            result[match_count++] = row_idx;
         }
      }

      return match_count;
   }

   uint32_t SkippingScanPerRowSimple(const FsstBlock& block, std::vector<uint32_t>& result, const std::vector<uint8_t>& required_symbols)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);
         bool has_any_of_the_required_symbol = std::any_of(required_symbols.begin(), required_symbols.end(), [&](uint8_t symbol) {
            return compressed_text.find(symbol) != std::string_view::npos;
         });
         if (!has_any_of_the_required_symbol) {
            skipped++;
            continue;
         }

         // Matching code.
         const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(compressed_text.data());
         bool match = state_machine.fsst_lookup_zerokmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         if (match) {
            result[match_count++] = row_idx;
         }
      }

      return match_count;
   }

   uint32_t SkippingScanPerRowSimd1(const FsstBlock& block, std::vector<uint32_t>& result, const std::vector<uint8_t>& required_symbols)
   {
      assert(required_symbols.size() == 1);
      uint8_t symbol = required_symbols[0];
      __m256i symbol_vec = _mm256_set1_epi8(symbol);

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Check for required symbol.
         uint32_t has_any_of_the_required_symbol = 0;
         std::string_view compressed_text = block.GetRow(row_idx);
         for (uint32_t idx = 0; idx < compressed_text.size(); idx += 32) {
            __m256i value_vec = _mm256_loadu_si256((__m256i*)(compressed_text.data() + idx));
            __m256i byte_mask_vec = _mm256_cmpeq_epi8(value_vec, symbol_vec);
            has_any_of_the_required_symbol |= _mm256_movemask_epi8(byte_mask_vec);
         }

         // No math -> easy win.
         if (!has_any_of_the_required_symbol) {
            skipped++;
            continue;
         }

         // Otherwise, use Mihail's slow but still fast matcher.
         const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(compressed_text.data());
         bool match = state_machine.fsst_lookup_zerokmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         if (match) {
            result[match_count++] = row_idx;
         }
      }

      return match_count;
   }

   uint32_t SkippingScanPerRowSimd2(const FsstBlock& block, std::vector<uint32_t>& result, const std::vector<uint8_t>& required_symbols)
   {
      assert(required_symbols.size() == 2);
      uint8_t symbol_0 = required_symbols[0];
      uint8_t symbol_1 = required_symbols[1];
      __m256i symbol_vec_0 = _mm256_set1_epi8(symbol_0);
      __m256i symbol_vec_1 = _mm256_set1_epi8(symbol_1);

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Check for required symbol.
         uint32_t has_any_of_the_required_symbol = 0;
         std::string_view compressed_text = block.GetRow(row_idx);
         for (uint32_t idx = 0; idx < compressed_text.size(); idx += 32) {
            __m256i value_vec = _mm256_loadu_si256((__m256i*)(block.data.data() + idx));
            __m256i byte_mask_vec = _mm256_or_si256(
                _mm256_cmpeq_epi8(value_vec, symbol_vec_0),
                _mm256_cmpeq_epi8(value_vec, symbol_vec_1));
            has_any_of_the_required_symbol |= _mm256_movemask_epi8(byte_mask_vec);
         }

         // No math -> easy win.
         if (!has_any_of_the_required_symbol) {
            skipped++;
            continue;
         }

         // Otherwise, use Mihail's slow but still fast matcher.
         const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(compressed_text.data());
         bool match = state_machine.fsst_lookup_zerokmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         if (match) {
            result[match_count++] = row_idx;
         }
      }

      return match_count;
   }

   std::vector<uint8_t> byte_mask;

   void CreateByteMask(const FsstBlock& block, const std::vector<uint8_t>& required_symbols)
   {
      // Reset byte mask.
      uint32_t text_size = block.data.size();
      byte_mask.resize(text_size);
      memset(byte_mask.data(), 0, text_size);

      // SIMD scan.
      if (required_symbols.size() == 1) {
         uint8_t symbol = required_symbols[0];
         __m256i symbol_vec = _mm256_set1_epi8(symbol);
         uint32_t idx = 0;
         for (; idx + 31 < text_size; idx += 32) {
            __m256i value_vec = _mm256_loadu_si256((__m256i*)(block.data.data() + idx));
            __m256i byte_mask_vec = _mm256_cmpeq_epi8(value_vec, symbol_vec);
            _mm256_storeu_si256((__m256i*)(byte_mask.data() + idx), byte_mask_vec);
         }
         for (; idx < text_size; idx++) {
            byte_mask[idx] = byte_mask[idx] == symbol ? 0xFF : 0;
         }
         return;
      }
      else if (required_symbols.size() == 2) {
         uint8_t symbol_0 = required_symbols[0];
         uint8_t symbol_1 = required_symbols[1];
         __m256i symbol_vec_0 = _mm256_set1_epi8(symbol_0);
         __m256i symbol_vec_1 = _mm256_set1_epi8(symbol_1);
         uint32_t idx = 0;
         for (; idx + 31 < text_size; idx += 32) {
            __m256i value_vec = _mm256_loadu_si256((__m256i*)(block.data.data() + idx));
            __m256i byte_mask_vec = _mm256_or_si256(
                _mm256_cmpeq_epi8(value_vec, symbol_vec_0),
                _mm256_cmpeq_epi8(value_vec, symbol_vec_1));
            _mm256_storeu_si256((__m256i*)(byte_mask.data() + idx), byte_mask_vec);
         }
         for (; idx < text_size; idx++) {
            byte_mask[idx] = (byte_mask[idx] == symbol_0 || byte_mask[idx] == symbol_1) ? 0xFF : 0;
         }
         return;
      }
      else {
         throw "More than 2 symbols in the byte mask not implemented.";
      }
   }

   uint32_t SkippingScanByteMask(const FsstBlock& block, std::vector<uint32_t>& result, const std::vector<uint8_t>& required_symbols)
   {
      CreateByteMask(block, required_symbols);

      uint32_t match_count = 0;
      uint32_t limit = block.data.size();
      uint32_t row_idx = 0;
      uint32_t pos = 0;
      while (row_idx < block.row_count) {
         // Seek to the next potential match.
         while (pos < limit && byte_mask[pos] == 0) {
            pos++;
         }

         // Find the row for this position.
         while (block.offsets[row_idx + 1] < pos) {
            row_idx++;
         }

         // Matching code.
         std::string_view compressed_text = block.GetRow(row_idx);
         const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(compressed_text.data());
         // bool match = state_machine.fsst_lookup_kmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         bool match = state_machine.fsst_lookup_zerokmp_match(block.decoder, compressed_text.size(), cast_input, block.decoder.GetIdealBufferSize(compressed_text.size()));
         if (match) {
            result[match_count++] = row_idx;
         }
         pos = block.offsets[row_idx];
         row_idx++;
      }

      return match_count;
   }

private:
   std::string_view pattern;
   std::vector<char> decode_buffer;
   StateMachine2 state_machine;

   void PrintPaths(const FsstBlock& block, const std::set<std::vector<uint8_t>>& possible_full_paths)
   {
      std::cout << "There are " << possible_full_paths.size() << " paths through the pattern." << std::endl;
      for (const auto& path : possible_full_paths) {
         std::cout << "path: ";
         for (uint8_t symbol : path) {
            std::cout << block.decoder.SymbolToStr(symbol) << " (" << (int)symbol << ") ";
         }
         std::cout << std::endl;
      }
   }

   // Entries from symbol table where their suffix is a prefix of the pattern.
   std::vector<uint8_t> GetSymbolsWithSuffix(const FsstBlock& block, uint32_t shift, std::string_view pattern) const
   {
      std::vector<uint8_t> symbols;
      for (uint32_t symbol = 0; symbol < block.decoder.GetSymbolTableSize(); symbol++) {
         std::string symbol_text = block.decoder.SymbolToStr(symbol);
         if (symbol_text.size() <= shift) {
            continue;
         }
         std::string prefix = symbol_text.substr(shift);
         if (!prefix.empty() && pattern.starts_with(prefix)) {
            symbols.push_back(symbol);
         }
      }

      return symbols;
   }

   std::vector<uint8_t> GetSymbolsWithPrefix(const FsstBlock& block, std::string_view pattern) const
   {
      std::vector<uint8_t> symbols;
      for (uint32_t symbol = 0; symbol < block.decoder.GetSymbolTableSize(); symbol++) {
         std::string text = block.decoder.SymbolToStr(symbol);
         if (text.starts_with(pattern) || pattern.starts_with(text)) {
            symbols.push_back(symbol);
         }
      }
      return symbols;
   }

   std::vector<uint8_t> GetSymbolsContainingFullPattern(const FsstBlock& block) const
   {
      std::vector<uint8_t> symbols;
      for (uint32_t symbol = 0; symbol < block.decoder.GetSymbolTableSize(); symbol++) {
         std::string text = block.decoder.SymbolToStr(symbol);
         if (text.find(pattern) != std::string::npos) {
            symbols.push_back(symbol);
         }
      }
      return symbols;
   }

   void BuildFullPaths(const FsstBlock& block, uint32_t idx, const std::vector<uint8_t>& path, std::set<std::vector<uint8_t>>& result)
   {
      auto is_possible = [&](std::string_view symbol_text, std::vector<uint8_t>& symbols, uint32_t remaining_pattern_size) {
         if (symbol_text.size() >= remaining_pattern_size) {
            return true;
         }
         for (auto symbol : symbols) {
            std::string other_symbol_text = block.decoder.SymbolToStr(symbol);
            if (symbol_text.size() < other_symbol_text.size() && other_symbol_text.size() <= remaining_pattern_size) {
               return false;
            }
         }
         return true;
      };

      std::vector<uint8_t> symbols = GetSymbolsWithPrefix(block, pattern.substr(idx));
      for (uint8_t symbol : symbols) {
         std::string text = block.decoder.SymbolToStr(symbol);
         bool possible = is_possible(text, symbols, pattern.size() - idx);
         if (possible) {
            std::vector<uint8_t> new_path = path;
            new_path.push_back(symbol);
            uint32_t new_idx = idx + text.size();
            if (new_idx >= pattern.size()) {
               result.insert(new_path);
            }
            else {
               BuildFullPaths(block, new_idx, new_path, result);
            }
         }
      }
   }

   std::vector<uint8_t> CreateRequiredSymbols(const FsstBlock& block)
   {
      std::set<std::vector<uint8_t>> possible_full_paths;

      // Symbols that contain the full pattern.
      {
         std::vector<uint8_t> symbols = GetSymbolsContainingFullPattern(block);
         for (uint8_t symbol : symbols) {
            std::string text = block.decoder.SymbolToStr(symbol);
            possible_full_paths.insert({symbol});
         }
      }

      // All symbols that have a suffix with which the pattern starts.
      uint32_t position_count = std::min(pattern.size(), size_t(8));
      std::vector<std::set<uint8_t>> symbols_to_get_to_starting_positions(position_count);
      for (uint32_t offset = 0; offset < 8; offset++) {
         std::vector<uint8_t> symbols = GetSymbolsWithSuffix(block, offset, pattern);
         for (uint8_t symbol : symbols) {
            std::string text = block.decoder.SymbolToStr(symbol);
            uint32_t starting_position = text.size() - offset;
            if (starting_position < symbols_to_get_to_starting_positions.size()) {
               symbols_to_get_to_starting_positions[starting_position].insert(symbol);
            }
         }
      }

      // auto is_possible = [&](std::string_view symbol_text, std::vector<uint8_t>& symbols, uint32_t remaining_pattern_size) {
      //    if (symbol_text.size() >= remaining_pattern_size) {
      //       return true;
      //    }
      //    for (auto symbol : symbols) {
      //       std::string text = block.decoder.SymbolToStr(symbol);
      //       if (text.size() > symbol_text.size()) {
      //          return false;
      //       }
      //    }
      //    return true;
      // };

      // Just print!!
      // std::cout << "--------------" << std::endl;
      // for (uint32_t idx = 0; idx < symbols_to_get_to_starting_positions.size(); idx++) {
      //    std::cout << "ways to get to " << idx << ": " << symbols_to_get_to_starting_positions[idx].size() << std::endl;
      //    if (symbols_to_get_to_starting_positions[idx].size() > 0 || idx == 0) {
      //       std::vector<uint8_t> symbols = GetSymbolsWithSuffix(block, 0, pattern.substr(idx));
      //       for (uint8_t symbol : symbols) {
      //          std::string text = block.decoder.SymbolToStr(symbol);
      //          bool possible = is_possible(text, symbols, pattern.size() - idx);
      //          std::cout << "symbol: " << (int)symbol << " text: " << text << " possible: " << possible << std::endl;
      //       }
      //    }
      // }

      for (uint32_t start_idx = 0; start_idx < symbols_to_get_to_starting_positions.size(); start_idx++) {
         if (start_idx == 0) {
            std::vector<uint8_t> path;
            BuildFullPaths(block, start_idx, path, possible_full_paths);
         }
         for (uint8_t symbol : symbols_to_get_to_starting_positions[start_idx]) {
            std::vector<uint8_t> path = {symbol};
            BuildFullPaths(block, start_idx, path, possible_full_paths);
         }
      }

      // Just print!!
      // PrintPaths(block, possible_full_paths);

      // Find the most common symbols.
      std::unordered_map<uint8_t, uint32_t> symbol_to_count;
      for (const auto& path : possible_full_paths) {
         for (uint8_t symbol : path) {
            symbol_to_count[symbol]++;
         }
      }
      std::vector<std::pair<uint8_t, uint32_t>> symbol_to_count_vec(symbol_to_count.begin(), symbol_to_count.end());
      std::sort(symbol_to_count_vec.begin(), symbol_to_count_vec.end(), [](const auto& a, const auto& b) {
         return a.second < b.second;
      });

      // Find a set of symbols that cover each path through the pattern.
      std::vector<uint8_t> required_symbols;
      while (!possible_full_paths.empty()) {
         uint8_t symbol = symbol_to_count_vec.back().first;
         symbol_to_count_vec.pop_back();
         bool seen = false;
         for (auto it = possible_full_paths.begin(); it != possible_full_paths.end();) {
            if (std::find(it->begin(), it->end(), symbol) != it->end()) {
               it = possible_full_paths.erase(it);
               seen = true;
            }
            else {
               ++it;
            }
         }
         if (seen) {
            required_symbols.push_back(symbol);
         }
      }

      return required_symbols;
   }
};
// -------------------------------------------------------------------------------------
class SkippingEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (std::count(pattern.begin(), pattern.end(), '%') == 2 &&
          pattern.find('_') == std::string::npos &&
          pattern.starts_with('%') &&
          pattern.ends_with('%')) {
         // Cut the pattern.
         auto cut_pattern = pattern.substr(1, pattern.size() - 2);

         // Build the state machine.
         StateMachine2 state_machine(cut_pattern);

         // Skip engine on top.
         return std::make_unique<SkippingEngine>(cut_pattern, std::move(state_machine));
      }
      return nullptr;
   }

   std::string GetName() final { return "skipping"; }
};
// -------------------------------------------------------------------------------------
