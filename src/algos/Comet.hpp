#pragma once
// -------------------------------------------------------------------------------------
#include "BenchmarkDriver.hpp"
#include "StateMachine.hpp"
#include <algorithm>
#include <span>
// -------------------------------------------------------------------------------------
template <class T>
class CometEngine : public Engine {
public:
   CometEngine(std::string_view pattern, StateMachine stateMachine)
   : pattern(pattern), stateMachine(std::move(stateMachine)) {}

   uint32_t Scan(const RawBlock& block, std::vector<uint32_t>& result)
   {
      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         if (static_cast<T*>(this)->RawMatches(block.GetRow(row_idx))) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result) final
   {
      // Init the FSST symbols.
      stateMachine.init_fsst_symbols(block.decoder);

      // Build the state lookup table.
      stateMachine.build_lookup_table();

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);

         // Match.
         if (!static_cast<T*>(this)->FsstMatches(block.decoder, compressed_text)) {
            // Decode the row.
            std::vector<char> decode_buffer;
            uint32_t ideal_buffer_size = block.decoder.GetIdealBufferSize(compressed_text.size());
            if (ideal_buffer_size > decode_buffer.size()) {
               decode_buffer.resize(ideal_buffer_size);
            }
            uint32_t decoded_size = block.decoder.Decode(compressed_text, decode_buffer);

            // Match.
            std::string_view text(decode_buffer.data(), decoded_size);
            if (text.find(pattern) != std::string_view::npos) {
               block.decoder.PrintSymbolTable(std::cerr);

               std::cerr << "row_idx=" << row_idx << std::endl;
               std::cerr << "text=" << text << std::endl;

               static_cast<T*>(this)->FsstMatches(block.decoder, compressed_text, true);

               stateMachine.init_fsst_symbols(block.decoder);
               assert(0);
            }
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

protected:
   std::string_view pattern;
   StateMachine stateMachine;
};
// -------------------------------------------------------------------------------------
class CometKmpEngine : public CometEngine<CometKmpEngine> {
public:
   CometKmpEngine(std::string_view pattern, StateMachine stateMachine)
       : CometEngine(pattern, std::move(stateMachine)) {}

   bool RawMatches(std::span<const char> input) noexcept {
      return stateMachine.kmp_match(input.data(), input.size());
   }

   bool FsstMatches(const FsstDecoder& fsstDecoder, std::span<const char> input, bool verbose = false) noexcept {
      const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(input.data());
      return stateMachine.fsst_kmp_match(fsstDecoder, input.size(), cast_input, fsstDecoder.GetIdealBufferSize(input.size()), verbose);
   }
};
// -------------------------------------------------------------------------------------
class CometEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      // Simple infix: %pattern%.
      if (std::count(pattern.begin(), pattern.end(), '%') == 2 && pattern.find('_') == std::string::npos &&
         pattern.starts_with('%') &&
         pattern.ends_with('%')) {
         // Cut the pattern.
         auto cut_pattern = pattern.substr(1, pattern.size() - 2);

         // Build the state machine.
         auto stateMachine = StateMachine(cut_pattern);
            
         // And return the engine.   
         return std::make_unique<CometKmpEngine>(cut_pattern, stateMachine);
      }

      return nullptr;
   }

   std::string GetName() final { return "comet[lookup-kmp]"; }
};
// -------------------------------------------------------------------------------------
