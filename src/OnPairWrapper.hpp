#pragma once
// -------------------------------------------------------------------------------------
#include <cstdint>
#include <string_view>
#include <vector>
#include <onpair/api.h>
#include "Utility.hpp"
// -------------------------------------------------------------------------------------
class OnPairColumnHandle : public NonCopyable {
public:
   OnPairColumnHandle() = default;
   OnPairColumnHandle(OnPairColumnHandle&&) = default;
   OnPairColumnHandle& operator=(OnPairColumnHandle&&) = default;

   void Build(const char* data, const uint32_t* offsets, std::size_t row_count,
              const onpair::encoding::TrainingConfig& cfg = {});

   onpair::OnPairColumnView view() const noexcept { return column_.view(); }
   std::size_t num_strings() const noexcept { return column_.num_strings(); }
   std::size_t bytes_used() const noexcept { return column_.bytes_used(); }

   std::uint32_t Decode(std::uint32_t row_idx, std::vector<char>& buf) const;

private:
   onpair::OnPairColumn column_;
};
// -------------------------------------------------------------------------------------
