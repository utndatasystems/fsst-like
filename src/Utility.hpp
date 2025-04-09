#pragma once
// -------------------------------------------------------------------------------------
#include <cassert>
#include <cstdint>
// -------------------------------------------------------------------------------------
struct NonCopyable {
   NonCopyable() = default;
   NonCopyable(const NonCopyable &) = delete;
   NonCopyable &operator=(const NonCopyable &) = delete;
   NonCopyable(NonCopyable &&) = default;
   NonCopyable &operator=(NonCopyable &&) = default;
};
// -------------------------------------------------------------------------------------
