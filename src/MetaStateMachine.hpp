#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <set>
#include "BenchmarkDriver.hpp"
#include "StateMachine.hpp"
#include "Utility.hpp"
// -------------------------------------------------------------------------------------
class MetaStateMachine {
public:
   MetaStateMachine(std::string_view pattern)
   {
      // TODO: Check if we need auto&& for `local_pattern`.
      auto local_patterns = SplitPattern(pattern);
      for (auto local_pattern : local_patterns) {
         state_machines.push_back(StateMachineView(local_pattern));
      }
      num_machines = state_machines.size();
   }

   void init(const FsstDecoder& fsstDecoder)
   {
      // TODO: Maybe reuse the fsst table, since it's anyway the same.
      for (unsigned index = 0; index != num_machines; ++index)
         state_machines[index].init(fsstDecoder);
   }

   void precompute() {
      for (unsigned index = 0; index != num_machines; ++index)
         state_machines[index].precompute();
   }

   bool fsst_lookup_kmp_match(const FsstDecoder& fsstDecoder, size_t lenIn, const unsigned char* strIn, size_t size)
   {
      unsigned machine_index = 0;
      state_machines[machine_index].init_state();

      // TODO: Remove the `verbose` part afterwards!!!
      auto consume_char = [this, &machine_index](unsigned char c) {
         // Accept the char.
         state_machines[machine_index].accept(c);

         // Not yet finished with this machine?
         if (state_machines[machine_index].get_state() != state_machines[machine_index].size())
            return true;

         // Otherwise, move to the next one.
         ++machine_index;

         // Done?
         if (machine_index == num_machines)
            return false;
         
         // Init the next machine.
         state_machines[machine_index].init_state();

         // We still need to go.
         return true;
      };

      auto consume_code = [this, &machine_index, consume_char](size_t code) {
         // TODO: Store the position where the machine accepted.
         // TODO: We need to continue from that part of the code.
         auto stop_pos = state_machines[machine_index].accept_symbol_with_lookup(code);

         // Not yet finished with this machine?
         if (state_machines[machine_index].get_state() != state_machines[machine_index].size())
            return true;

         // Otherwise, move to the next one.
         ++machine_index;

         // Was this the last one?
         if (machine_index == num_machines)
            return false;
         
         // Init the next machine.
         state_machines[machine_index].init_state();

         // Note: This is still very naive.
         auto& fsst_symbols = state_machines[machine_index].fsst_symbols;
         for (unsigned index = stop_pos, limit = fsst_symbols[code].size(); index != limit; ++index) {
            // Note: We _do_ have to use `consume_char`.
            // Note: But you can pre-compute how much from the current symbol the current machine "eats".
            if (!consume_char(fsst_symbols[code][index]))
               return false;
         }

         // We still need to go.
         return true;
      };

      return fsstDecoder.Iterate(lenIn, strIn, size, consume_char, consume_code);
   }

private:
   unsigned num_machines;
   std::vector<StateMachineView> state_machines;
};