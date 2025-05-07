#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <set>
#include "BenchmarkDriver.hpp"
#include "Comet.hpp"
#include "Utility.hpp"
// -------------------------------------------------------------------------------------
class SkippingEngine : public Engine {
public:
   SkippingEngine(std::string_view pattern, std::unique_ptr<CometEngine<CometKmpEngine>> comet)
       : pattern(pattern)
       , decode_buffer(128)
       , comet(move(comet))
   {
   }

   ~SkippingEngine() override
   {
      for (const auto& [key, value] : required_symbol_count_agg) {
         std::cout << "required symbol count: " << key << " count: " << value << std::endl;
      }
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

   void BuildCorePaths(const FsstBlock& block, uint32_t idx, const std::vector<uint8_t>& path, std::set<std::vector<uint8_t>>& result)
   {
      auto is_possible = [&](std::string_view symbol_text, std::vector<uint8_t>& symbols, uint32_t remaining_pattern_size) {
         if (symbol_text.size() >= remaining_pattern_size) {
            return true;
         }
         for (auto symbol : symbols) {
            std::string text = block.decoder.SymbolToStr(symbol);
            if (text.size() > symbol_text.size()) {
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
            uint32_t new_idx = idx + text.size();
            if (new_idx >= pattern.size()) {
               if (new_path.empty()) {
                  new_path.push_back(symbol);
               }
               result.insert(new_path);
            }
            else {
               new_path.push_back(symbol);
               BuildCorePaths(block, new_idx, new_path, result);
            }
         }
      }
   }

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result) final
   {
      comet->InitializeForCompressedScan(block);

      // if (counter != 16) {
      //    counter++;
      //    return WordlessScan(block, result);
      // }

      std::cout << "=======================================================================" << std::endl;
      std::cout << "block_idx: " << counter << std::endl;
      counter++;
      block.decoder.PrintSymbolTable(std::cout);
      block.PrintUsedChars(std::cout);
      std::set<std::vector<uint8_t>> possible_full_paths;
      std::set<std::vector<uint8_t>> possible_paths;

      // Symbols that contain the full pattern.
      {
         std::cout << "complete matches:" << std::endl;
         std::vector<uint8_t> symbols = GetSymbolsContainingFullPattern(block);
         for (uint8_t symbol : symbols) {
            std::string text = block.decoder.SymbolToStr(symbol);
            std::cout << "symbol: " << (int)symbol << " text: " << text << std::endl;
            possible_paths.insert({symbol});
            possible_full_paths.insert({symbol});
         }
      }

      // All symbols that have a suffix with which the pattern starts.
      std::cout << "--------------" << std::endl;

      uint32_t position_count = std::min(pattern.size(), size_t(8));
      std::vector<std::set<uint8_t>> symbols_to_get_to_starting_positions(position_count);
      for (uint32_t offset = 0; offset < 8; offset++) {
         std::cout << "offset=" << offset << std::endl;
         std::vector<uint8_t> symbols = GetSymbolsWithSuffix(block, offset, pattern);
         for (uint8_t symbol : symbols) {
            std::string text = block.decoder.SymbolToStr(symbol);
            uint32_t starting_position = text.size() - offset;
            if (starting_position < symbols_to_get_to_starting_positions.size()) {
               symbols_to_get_to_starting_positions[starting_position].insert(symbol);
               std::cout << "symbol: " << (int)symbol << " text: " << text << std::endl;
            }
         }
      }

      auto is_possible = [&](std::string_view symbol_text, std::vector<uint8_t>& symbols, uint32_t remaining_pattern_size) {
         if (symbol_text.size() >= remaining_pattern_size) {
            return true;
         }
         for (auto symbol : symbols) {
            std::string text = block.decoder.SymbolToStr(symbol);
            if (text.size() > symbol_text.size()) {
               return false;
            }
         }
         return true;
      };

      // Just print!!
      std::cout << "--------------" << std::endl;
      for (uint32_t idx = 0; idx < symbols_to_get_to_starting_positions.size(); idx++) {
         std::cout << "ways to get to " << idx << ": " << symbols_to_get_to_starting_positions[idx].size() << std::endl;
         if (symbols_to_get_to_starting_positions[idx].size() > 0 || idx == 0) {
            std::vector<uint8_t> symbols = GetSymbolsWithSuffix(block, 0, pattern.substr(idx));
            for (uint8_t symbol : symbols) {
               std::string text = block.decoder.SymbolToStr(symbol);
               bool possible = is_possible(text, symbols, pattern.size() - idx);
               std::cout << "symbol: " << (int)symbol << " text: " << text << " possible: " << possible << std::endl;
            }
         }
      }

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
      std::cout << "There are " << possible_full_paths.size() << " paths through the pattern." << std::endl;
      for (const auto& path : possible_full_paths) {
         std::cout << "path: ";
         for (uint8_t symbol : path) {
            std::cout << block.decoder.SymbolToStr(symbol) << " (" << (int)symbol << ") ";
         }
         std::cout << std::endl;
      }

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

      // Print
      std::cout << "Required symbol count: " << required_symbols.size() << std::endl;
      required_symbol_count_agg[required_symbols.size()] += 1;
      for (uint8_t symbol : required_symbols) {
         std::cout << "required symbol: " << (int)symbol << " text: " << block.decoder.SymbolToStr(symbol) << std::endl;
      }

      // for (uint32_t start_idx = 0; start_idx < symbols_to_get_to_starting_positions.size(); start_idx++) {
      //    if (start_idx == 0) {
      //       std::vector<uint8_t> path;
      //       BuildCorePaths(block, start_idx, path, possible_paths);
      //    }
      //    if (!symbols_to_get_to_starting_positions[start_idx].empty()) {
      //       std::vector<uint8_t> path;
      //       BuildCorePaths(block, start_idx, path, possible_paths);
      //    }
      // }

      // // Just print!!
      // std::cout << "There are " << possible_paths.size() << " core paths through the pattern." << std::endl;
      // for (const auto& path : possible_paths) {
      //    std::cout << "path: ";
      //    for (uint8_t symbol : path) {
      //       std::cout << (int)symbol << " ";
      //    }
      //    std::cout << std::endl;
      // }

      std::cout << "scanning ... " << std::endl;

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);
         bool has_any_of_the_required_symbol = std::any_of(required_symbols.begin(), required_symbols.end(), [&](uint8_t symbol) {
            return compressed_text.find(symbol) != std::string_view::npos;
         });
         if (!has_any_of_the_required_symbol) {
            continue;
         }

         // Match.
         if (static_cast<CometKmpEngine*>(comet.get())->FsstMatches(block.decoder, compressed_text)) {
            result[match_count++] = row_idx;
            // if (!has_any_of_the_required_symbol) {
            //    std::cout << "row_idx: " << row_idx << " compressed_text: ";
            //    for (uint8_t symbol : compressed_text) {
            //       std::cout << (int)symbol << " ";
            //    }
            //    std::cout << "text: " << text << std::endl;
            //    std::cout << std::endl;
            // }
         }
      }

      std::cout << "all good!" << std::endl;
      return match_count;
   }

   std::unordered_map<uint32_t, uint32_t> required_symbol_count_agg;

   uint32_t WordlessScan(const FsstBlock& block, std::vector<uint32_t>& result)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);

         // Decode the row.
         uint32_t ideal_buffer_size = block.decoder.GetIdealBufferSize(compressed_text.size());
         if (ideal_buffer_size > decode_buffer.size()) {
            decode_buffer.resize(ideal_buffer_size);
         }
         uint32_t decoded_size = block.decoder.Decode(compressed_text, decode_buffer);

         // Match.
         std::string_view text(decode_buffer.data(), decoded_size);
         if (text.find(pattern) != std::string_view::npos) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

private:
   std::string_view pattern;
   std::vector<char> decode_buffer;
   std::unique_ptr<CometEngine<CometKmpEngine>> comet;

   int counter = 0;
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
         auto stateMachine = StateMachine(cut_pattern);

         // And return the engine.
         auto comet = std::make_unique<CometKmpEngine>(cut_pattern, stateMachine);

         // Skip engine on top.
         return std::make_unique<SkippingEngine>(cut_pattern, move(comet));
      }
      return nullptr;
   }

   std::string GetName() final { return "skipping"; }
};
// -------------------------------------------------------------------------------------
