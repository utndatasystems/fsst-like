#pragma once
// -------------------------------------------------------------------------------------
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <onpair/api.h>
#include "BenchmarkDriver.hpp"
#include "Utility.hpp"
// -------------------------------------------------------------------------------------
class OnPairContainsEngine : public Engine {
public:
   explicit OnPairContainsEngine(std::string_view p) : pattern(p) {}

   uint32_t Scan(const RawBlock& blk, std::vector<uint32_t>& res) override
   {
      uint32_t n = 0;
      for (uint32_t i = 0; i < blk.row_count; ++i) {
         if (blk.GetRow(i).find(pattern) != std::string_view::npos) res[n++] = i;
      }
      return n;
   }

   uint32_t Scan(const OnPairBlock& blk, std::vector<uint32_t>& res) override
   {
      auto v = blk.column.view();
      uint32_t n = 0;
      v.contains(pattern, [&](std::size_t idx) { res[n++] = static_cast<uint32_t>(idx); });
      return n;
   }

private:
   std::string pattern;
};
// -------------------------------------------------------------------------------------
class OnPairStartsWithEngine : public Engine {
public:
   explicit OnPairStartsWithEngine(std::string_view p) : pattern(p) {}

   uint32_t Scan(const RawBlock& blk, std::vector<uint32_t>& res) override
   {
      uint32_t n = 0;
      for (uint32_t i = 0; i < blk.row_count; ++i) {
         if (blk.GetRow(i).starts_with(pattern)) res[n++] = i;
      }
      return n;
   }

   uint32_t Scan(const OnPairBlock& blk, std::vector<uint32_t>& res) override
   {
      auto v = blk.column.view();
      uint32_t n = 0;
      v.starts_with(pattern, [&](std::size_t idx) { res[n++] = static_cast<uint32_t>(idx); });
      return n;
   }

private:
   std::string pattern;
};
// -------------------------------------------------------------------------------------
class OnPairEqualsEngine : public Engine {
public:
   explicit OnPairEqualsEngine(std::string_view p) : pattern(p) {}

   uint32_t Scan(const RawBlock& blk, std::vector<uint32_t>& res) override
   {
      uint32_t n = 0;
      for (uint32_t i = 0; i < blk.row_count; ++i) {
         if (blk.GetRow(i) == pattern) res[n++] = i;
      }
      return n;
   }

   uint32_t Scan(const OnPairBlock& blk, std::vector<uint32_t>& res) override
   {
      auto v = blk.column.view();
      uint32_t n = 0;
      v.equals(pattern, [&](std::size_t idx) { res[n++] = static_cast<uint32_t>(idx); });
      return n;
   }

private:
   std::string pattern;
};
// -------------------------------------------------------------------------------------
class OnPairMultiContainsEngine : public Engine {
public:
   explicit OnPairMultiContainsEngine(std::vector<std::string> ps)
       : patterns(std::move(ps)) {}

   uint32_t Scan(const RawBlock& blk, std::vector<uint32_t>& res) override
   {
      uint32_t n = 0;
      for (uint32_t i = 0; i < blk.row_count; ++i) {
         auto t = blk.GetRow(i);
         std::size_t pos = 0;
         bool ok = true;
         for (auto& p : patterns) {
            auto f = t.find(p, pos);
            if (f == std::string_view::npos) { ok = false; break; }
            pos = f + p.size();
         }
         if (ok) res[n++] = i;
      }
      return n;
   }

   uint32_t Scan(const OnPairBlock& blk, std::vector<uint32_t>& res) override
   {
      // OnPair's Aho-Corasick automaton has "any pattern matches" semantics; the
      // LIKE multi-substring '%p1%p2%...%' requires *ordered sequential* matches.
      // Decompress-then-match keeps this correct; an automaton-only path would
      // need conjoined per-pattern state tracking with position constraints.
      auto v = blk.column.view();
      std::vector<char> buf(4096 + onpair::DECOMPRESS_BUFFER_PADDING);
      uint32_t n = 0;
      const std::size_t total = v.num_strings();
      for (std::size_t i = 0; i < total; ++i) {
         auto len = v.decompress(i, buf.data());
         std::string_view t(buf.data(), len);
         std::size_t pos = 0;
         bool ok = true;
         for (auto& p : patterns) {
            auto f = t.find(p, pos);
            if (f == std::string_view::npos) { ok = false; break; }
            pos = f + p.size();
         }
         if (ok) res[n++] = static_cast<uint32_t>(i);
      }
      return n;
   }

private:
   std::vector<std::string> patterns;
};
// -------------------------------------------------------------------------------------
// Bulk-decompress + std baseline. Mirrors FSST's `stl` engine: decompress every
// row of the block once via OnPair's decompress_all, then run std::find /
// starts_with / ends_with over the flat buffer. Useful as an apples-to-apples
// comparison against FSST's `stl` and to give OnPair a suffix path.
template <class T>
class OnPairStdEngine : public Engine {
public:
   explicit OnPairStdEngine(std::string_view p) : pattern(p) {}

   uint32_t Scan(const RawBlock& blk, std::vector<uint32_t>& res) override
   {
      uint32_t n = 0;
      for (uint32_t i = 0; i < blk.row_count; ++i) {
         if (static_cast<T*>(this)->Matches(blk.GetRow(i))) res[n++] = i;
      }
      return n;
   }

   uint32_t Scan(const OnPairBlock& blk, std::vector<uint32_t>& res) override
   {
      auto v = blk.column.view();
      const std::size_t need = std::size_t{blk.decompressed_bytes} + onpair::DECOMPRESS_BUFFER_PADDING;
      if (decode_buffer.size() < need) decode_buffer.resize(need);
      if (decoded_offsets.size() < std::size_t{blk.row_count} + 1) decoded_offsets.resize(std::size_t{blk.row_count} + 1);

      v.decompress_all(decode_buffer.data(), decoded_offsets.data());

      uint32_t n = 0;
      for (uint32_t i = 0; i < blk.row_count; ++i) {
         std::string_view t(decode_buffer.data() + decoded_offsets[i],
                            decoded_offsets[i + 1] - decoded_offsets[i]);
         if (static_cast<T*>(this)->Matches(t)) res[n++] = i;
      }
      return n;
   }

protected:
   std::string pattern;
   std::vector<char> decode_buffer;
   std::vector<uint32_t> decoded_offsets;
};
// -------------------------------------------------------------------------------------
class OnPairStdContainsEngine : public OnPairStdEngine<OnPairStdContainsEngine> {
public:
   using OnPairStdEngine::OnPairStdEngine;
   bool Matches(std::string_view t) const noexcept { return t.find(pattern) != std::string_view::npos; }
};
// -------------------------------------------------------------------------------------
class OnPairStdStartsWithEngine : public OnPairStdEngine<OnPairStdStartsWithEngine> {
public:
   using OnPairStdEngine::OnPairStdEngine;
   bool Matches(std::string_view t) const noexcept { return t.starts_with(pattern); }
};
// -------------------------------------------------------------------------------------
class OnPairStdEndsWithEngine : public OnPairStdEngine<OnPairStdEndsWithEngine> {
public:
   using OnPairStdEngine::OnPairStdEngine;
   bool Matches(std::string_view t) const noexcept { return t.ends_with(pattern); }
};
// -------------------------------------------------------------------------------------
class OnPairStdGeneralEngine : public OnPairStdEngine<OnPairStdGeneralEngine> {
public:
   OnPairStdGeneralEngine(std::vector<std::string> ps)
       : OnPairStdEngine<OnPairStdGeneralEngine>({}), patterns_(std::move(ps)) {}

   bool Matches(std::string_view t) const noexcept
   {
      std::size_t pos = 0;
      for (auto& p : patterns_) {
         auto f = t.find(p, pos);
         if (f == std::string_view::npos) return false;
         pos = f + p.size();
      }
      return true;
   }

private:
   std::vector<std::string> patterns_;
};
// -------------------------------------------------------------------------------------
class OnPairStdEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (pattern.find('_') != std::string_view::npos) return nullptr;

      const bool sw = pattern.starts_with('%');
      const bool ew = pattern.ends_with('%');
      const auto pcount = std::count(pattern.begin(), pattern.end(), '%');

      if (pcount == 2 && sw && ew) {
         return std::make_unique<OnPairStdContainsEngine>(pattern.substr(1, pattern.size() - 2));
      }
      if (pcount == 1 && ew && !sw) {
         return std::make_unique<OnPairStdStartsWithEngine>(pattern.substr(0, pattern.size() - 1));
      }
      if (pcount == 1 && sw && !ew) {
         return std::make_unique<OnPairStdEndsWithEngine>(pattern.substr(1));
      }
      if (sw && ew) {
         auto inner = pattern.substr(1, pattern.size() - 2);
         auto views = SplitPattern(inner);
         std::vector<std::string> ps;
         ps.reserve(views.size());
         for (auto v : views) ps.emplace_back(v);
         return std::make_unique<OnPairStdGeneralEngine>(std::move(ps));
      }
      return nullptr;
   }

   std::string GetName() final { return "onpair-stl"; }
};
// -------------------------------------------------------------------------------------
class OnPairUnsupportedEngine : public Engine {
public:
   explicit OnPairUnsupportedEngine(std::string reason) : reason_(std::move(reason)) {}

   uint32_t Scan(const RawBlock&, std::vector<uint32_t>&) override
   {
      throw std::logic_error(reason_);
   }
   uint32_t Scan(const FsstBlock&, std::vector<uint32_t>&) override
   {
      throw std::logic_error(reason_);
   }
   uint32_t Scan(const OnPairBlock&, std::vector<uint32_t>&) override
   {
      throw std::logic_error(reason_);
   }

private:
   std::string reason_;
};
// -------------------------------------------------------------------------------------
class OnPairEngineFactory : public EngineFactory {
public:
   std::unique_ptr<Engine> Create(std::string_view pattern) final
   {
      if (pattern.find('_') != std::string_view::npos) return nullptr;

      const bool sw = pattern.starts_with('%');
      const bool ew = pattern.ends_with('%');
      const auto pcount = std::count(pattern.begin(), pattern.end(), '%');

      // No wildcards => equality.
      if (pcount == 0) {
         return std::make_unique<OnPairEqualsEngine>(pattern);
      }

      // 'p%' -> prefix
      if (pcount == 1 && ew && !sw) {
         return std::make_unique<OnPairStartsWithEngine>(pattern.substr(0, pattern.size() - 1));
      }

      // '%p' -> suffix => unsupported by OnPair.
      if (pcount == 1 && sw && !ew) {
         return std::make_unique<OnPairUnsupportedEngine>(
             "OnPair: suffix matching ('%p') is not implemented");
      }

      // '%p%' -> single substring
      if (pcount == 2 && sw && ew) {
         return std::make_unique<OnPairContainsEngine>(
             pattern.substr(1, pattern.size() - 2));
      }

      // '%p1%p2%...%' -> multi-substring (sequential)
      if (sw && ew) {
         auto inner = pattern.substr(1, pattern.size() - 2);
         auto views = SplitPattern(inner);
         std::vector<std::string> ps;
         ps.reserve(views.size());
         for (auto v : views) ps.emplace_back(v);
         return std::make_unique<OnPairMultiContainsEngine>(std::move(ps));
      }

      return nullptr;
   }

   std::string GetName() final { return "onpair"; }
};
// -------------------------------------------------------------------------------------
