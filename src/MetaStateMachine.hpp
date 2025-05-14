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
      auto local_patterns = SplitPattern(pattern);
      for (auto local_pattern : local_patterns) {
         state_machines.push_back(StateMachineView(local_pattern));
      }
      num_machines = state_machines.size();
   }

   void init()

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

      auto consume_code = [this, &machine_index](size_t code) {
         // TODO: Store the position where the machine accepted.
         // TODO: We need to continue from that part of the code.
         state_machines[machine_index].accept_symbol_with_lookup(code);

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

      return fsstDecoder.Iterate(lenIn, strIn, size, consume_char, consume_code);
   }

private:
   unsigned num_machines;
   std::vector<StateMachineView> state_machines;
};