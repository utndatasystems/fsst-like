#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <span>
#include "MetaStateMachine.hpp"
#include "BenchmarkDriver.hpp"
#include "StateMachine.hpp"
// -------------------------------------------------------------------------------------
template <class T, class MachineType>
class CometEngine : public Engine {
public:
   CometEngine(std::string_view pattern, MachineType machine)
       : pattern(pattern), machine(std::move(machine)) {}

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

   void InitializeForCompressedScan(const FsstBlock& block)
   {
      // Init the FSST symbols.
      machine.init(block.decoder);

      // Precompute the state lookup tables.
      machine.precompute();
   }

   uint32_t Scan(const FsstBlock& block, std::vector<uint32_t>& result) final
   {
      InitializeForCompressedScan(block);

      uint32_t match_count = 0;
      for (uint32_t row_idx = 0; row_idx < block.row_count; row_idx++) {
         // Get encoded row.
         std::string_view compressed_text = block.GetRow(row_idx);

         // Match.
         if (static_cast<T*>(this)->FsstMatches(block.decoder, compressed_text)) {
            result[match_count++] = row_idx;
         }
      }
      return match_count;
   }

protected:
   std::string_view pattern;
   MachineType machine;
};
// -------------------------------------------------------------------------------------
class CometKmpEngine : public CometEngine<CometKmpEngine, StateMachine> {
public:
   CometKmpEngine(std::string_view pattern, StateMachine machine)
       : CometEngine(pattern, std::move(machine)) {}

   bool RawMatches(std::string_view input) noexcept
   {
      return machine.raw_kmp_match(input.data(), input.size());
   }

   bool FsstMatches(const FsstDecoder& fsstDecoder, std::string_view input) noexcept
   {
      const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(input.data());
      return machine.fsst_lookup_kmp_match(fsstDecoder, input.size(), cast_input, fsstDecoder.GetIdealBufferSize(input.size()));
   }
};
// -------------------------------------------------------------------------------------
class CometKmpMetaEngine : public CometEngine<CometKmpMetaEngine, MetaStateMachine> {
public:
   CometKmpMetaEngine(std::string_view pattern, MetaStateMachine machine)
       : CometEngine(pattern, std::move(machine)) {}

   bool RawMatches(std::string_view input) noexcept
   {
      // return stateMachine.kmp_match(input.data(), input.size());
   }

   bool FsstMatches(const FsstDecoder& fsstDecoder, std::string_view input) noexcept
   {
      const unsigned char* cast_input = reinterpret_cast<const unsigned char*>(input.data());
      return machine.fsst_lookup_kmp_match(fsstDecoder, input.size(), cast_input, fsstDecoder.GetIdealBufferSize(input.size()));
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

         // // Can we enable ZeroKMP?
         // if (stateMachine.has_zerokmp_property()) {
         //    return std::make_unique<CometSingleZeroKmpEngine>(cut_pattern, stateMachine);
         // }

         // And return the engine.
         return std::make_unique<CometKmpEngine>(cut_pattern, stateMachine);
      }

      // Generic case: %pattern1%pattern2%.
      // TODO: Implement: [%]pattern1%pattern2[%].
      if (pattern.starts_with('%') &&
         pattern.starts_with('%') &&
         pattern.find('_') == std::string::npos) {
         auto meta_state_machine = MetaStateMachine(pattern);

         return std::make_unique<CometKmpMetaEngine>(pattern, meta_state_machine);
      }

      return nullptr;
   }

   std::string GetName() final { return "comet[lookup-kmp]"; }
};
// -------------------------------------------------------------------------------------
