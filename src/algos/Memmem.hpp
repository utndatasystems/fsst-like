#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <ranges>
#include "BenchmarkDriver.hpp"
#include "StdFind.hpp"
// -------------------------------------------------------------------------------------
class GeneralMemmemEngine : public StdEngine<GeneralMemmemEngine> {
public:
    explicit GeneralMemmemEngine(std::string_view pattern)
        : StdEngine(pattern) 
    {
        patterns_ = SplitPattern(pattern);
    }

    bool Matches(std::string_view text) const noexcept {
        std::size_t start_pos = 0;
        for (const auto& pat : patterns_) {
            if (text.size() <= start_pos)
                return false;

            const void* match = memmem(
                text.data() + start_pos,
                text.size() - start_pos,
                pat.data(),
                pat.size()
            );

            if (!match)
                return false;

            // Advance `start_pos` to after the matched pattern.
            start_pos = static_cast<const char*>(match) - text.data() + pat.size();
        }
        return true;
    }

private:
    std::vector<std::string_view> patterns_;
};
// -------------------------------------------------------------------------------------
class MemmemEngine : public StdEngine<MemmemEngine> {
public:
   MemmemEngine(std::string_view pattern)
       : StdEngine(pattern) {}

   bool Matches(std::string_view text) noexcept {
        const void* match = memmem(
            text.data(),
            text.size(),
            pattern.data(),
            pattern.size()
        );

        return (!!match);
    }
};
// -------------------------------------------------------------------------------------
class MemmemEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (std::count(pattern.begin(), pattern.end(), '%') == 2 &&
          pattern.find('_') == std::string::npos &&
          pattern.starts_with('%') &&
          pattern.ends_with('%')) {
         return std::make_unique<MemmemEngine>(pattern.substr(1, pattern.size() - 2));
      }
      if (std::count(pattern.begin(), pattern.end(), '%') == 1 &&
          pattern.find('_') == std::string::npos &&
          pattern.ends_with('%')) {
         return std::make_unique<StdStartsWithEngine>(pattern.substr(0, pattern.size() - 1));
      }
      if (std::count(pattern.begin(), pattern.end(), '%') == 1 &&
          pattern.find('_') == std::string::npos &&
          pattern.starts_with('%')) {
         return std::make_unique<StdEndsWithEngine>(pattern.substr(1, pattern.size() - 1));
      }

      // Generic case for: %p1%p2%.
      // TODO: We also need the more general case: [%]p1%p2[%].
      if (pattern.starts_with('%') &&
          pattern.ends_with('%') &&
          pattern.find('_') == std::string::npos) {
         return std::make_unique<GeneralMemmemEngine>(pattern.substr(1, pattern.size() - 2));
      }

      return nullptr;
   }

   std::string GetName() final { return "memmem"; }
};
// -------------------------------------------------------------------------------------
