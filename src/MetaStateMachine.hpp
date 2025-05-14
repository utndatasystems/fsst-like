#pragma once
// -------------------------------------------------------------------------------------
#include "StateMachine.hpp"
// -------------------------------------------------------------------------------------
static constexpr bool COMPLETED = false;
// -------------------------------------------------------------------------------------
class MetaStateMachine {
public:
   MetaStateMachine(std::string_view pattern)
   {
      // TODO: Check if we need auto&& for `local_pattern`.
      auto local_patterns = SplitPattern(pattern);
      for (auto local_pattern : local_patterns) {
         std::cerr << "local_pattern=" << local_pattern << std::endl;
         state_machines.push_back(StateMachineView(local_pattern));
      }
      num_machines = state_machines.size();
   }

   void init(const FsstDecoder& fsstDecoder)
   {
      std::cerr << "[init] start" << std::endl;
      // TODO: Maybe reuse the fsst table, since it's anyway the same.
      for (unsigned index = 0; index != num_machines; ++index)
         state_machines[index].init(fsstDecoder);
      std::cerr << "[init] stop" << std::endl;
   }

   void precompute() {
      std::cerr << "[precompute] start" << std::endl;
      for (unsigned index = 0; index != num_machines; ++index)
         state_machines[index].precompute();
      std::cerr << "[precompute] stop" << std::endl;
   }

   // TODO: Remove `verbose`.
   bool fsst_lookup_kmp_match(const FsstDecoder& fsstDecoder, size_t lenIn, const unsigned char* strIn, size_t size, bool verbose = false)
   {
      unsigned machine_index = 0;
      state_machines[machine_index].init_state();

      // TODO: Remove the `verbose` part afterwards!!!
      auto consume_char = [this, &machine_index, verbose](unsigned char c) {
         if (verbose) {
            std::cerr << "[consume_char] c=" << c << std::endl;
         }

         // Accept the char.
         state_machines[machine_index].accept(c);

         // Not yet finished with this machine?
         if (state_machines[machine_index].get_state() != state_machines[machine_index].size())
            return !COMPLETED;

         // Otherwise, move to the next one.
         ++machine_index;

         // Done?
         if (machine_index == num_machines)
            return COMPLETED;
         
         // Init the next machine.
         state_machines[machine_index].init_state();

         // We still need to go.
         return !COMPLETED;
      };

      auto consume_code = [this, &machine_index, consume_char, verbose](size_t code) {
         if (verbose) {
            std::cerr << "[consume_code] code=" << code << std::endl;
         }
         // TODO: Store the position where the machine accepted.
         // TODO: We need to continue from that part of the code.
         auto stop_pos = state_machines[machine_index].accept_symbol_with_lookup(code);

         if (verbose) {
            std::cerr << "stop_pos=" << stop_pos << std::endl;
         }

         // Not yet finished with this machine?
         if (state_machines[machine_index].get_state() != state_machines[machine_index].size())
            return !COMPLETED;

         // Otherwise, move to the next one.
         ++machine_index;

         // Was this the last one?
         if (machine_index == num_machines)
            return COMPLETED;
         
         // Init the next machine.
         state_machines[machine_index].init_state();

         // Note: This is still very naive.
         auto& fsst_symbols = state_machines[machine_index].fsst_symbols;
         for (unsigned index = stop_pos, limit = fsst_symbols[code].size(); index != limit; ++index) {
            // Note: We _do_ have to use `consume_char`.
            // Note: But you can pre-compute how much from the current symbol the current machine "eats".
            if (consume_char(fsst_symbols[code][index]) == COMPLETED)
               return COMPLETED;
         }

         // We still need to go.
         return !COMPLETED;
      };

      return fsstDecoder.Iterate(lenIn, strIn, size, consume_char, consume_code);
   }

private:
   unsigned num_machines;
   std::vector<StateMachineView> state_machines;
};