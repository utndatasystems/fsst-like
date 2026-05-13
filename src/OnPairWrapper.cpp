#include "OnPairWrapper.hpp"
// -------------------------------------------------------------------------------------
void OnPairColumnHandle::Build(const char* data, const uint32_t* offsets,
                               std::size_t row_count,
                               const onpair::encoding::TrainingConfig& cfg)
{
   column_ = onpair::OnPairColumn::compress(data, offsets, row_count, cfg);
}
// -------------------------------------------------------------------------------------
std::uint32_t OnPairColumnHandle::Decode(std::uint32_t row_idx,
                                         std::vector<char>& buf) const
{
   auto v = column_.view();
   const std::size_t need = 4096 + onpair::DECOMPRESS_BUFFER_PADDING;
   if (buf.size() < need) buf.resize(need);
   return static_cast<std::uint32_t>(v.decompress(row_idx, buf.data()));
}
// -------------------------------------------------------------------------------------
