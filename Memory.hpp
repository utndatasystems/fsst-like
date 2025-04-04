#pragma once
// -------------------------------------------------------------------------------------
#include <array>
#include <cstring>
#include <cstdint>
// -------------------------------------------------------------------------------------
struct Memory {
   // Store bytes at the given pointer.
   template <class T>
   static void Store(uint8_t *data, T value)
   {
      memcpy(data, &value, sizeof(T));
   }

   // Unaligned 4 byte load (assumes `data` is long enough).
   template <class RESULT, class T>
   static RESULT Load(const T *data)
   {
      static_assert(std::is_same_v<T, uint8_t> || std::is_same_v<T, char>);
      RESULT result;
      memcpy(&result, data, sizeof(RESULT));
      return result;
   }

   static constexpr uint64_t Broadcast(uint8_t num) { return uint64_t(num) << 56u | uint64_t(num) << 48u | uint64_t(num) << 40u | uint64_t(num) << 32u | uint64_t(num) << 24u | uint64_t(num) << 16u | uint64_t(num) << 8u | uint64_t(num) << 0u; }

   // Unaligned 8 byte load from `data` with size limit.
   _Pragma("GCC diagnostic push")
   _Pragma("GCC diagnostic ignored \"-Warray-bounds\"")
   template <class T>
   static uint64_t __attribute__((no_sanitize_address, no_sanitize_undefined)) LoadSafe8(const T *data, uint32_t size)
   {
      static_assert(std::is_same_v<T, uint8_t> || std::is_same_v<T, char>);

      // Shortcuts.
      if (size == 0) {
         return 0;
      }
      if (size >= 8) {
         return Load<uint64_t>(data);
      }

#ifdef VALGRIND_BUILD
      {
         // Valgrind is so not chill about this hacky code -> just memcpy instead.
         uint64_t result = 0;
         memcpy(&result, data, size);
         return result;
      }
#endif

      // We want to load 8 bytes: figure out where on the page `data` starts.
      uint64_t ptr = reinterpret_cast<uint64_t>(data);
      uint64_t cache_line_offset = ptr & 4095;

      // For the shifting below we assume little endian.

      // If `data` starts early enough on the page to contain `data+8` -> safe to load directly.
      if (cache_line_offset <= 4090) {
         uint64_t block = DirtyLoad(data);
         uint64_t mask = (1ull << (size * 8)) - 1;
         return block & mask;
      }

      // Otherwise, the string starts towards the end of the page.
      // Hence, it _might_ leak into the next one, depending on `size`.
      // Make sure we only touch the next page, if `data` actually spans there.
      uint32_t offset = (8 - size);
      uint64_t block = DirtyLoad(data - offset);
      return block >> (offset * 8);
   }
   _Pragma("GCC diagnostic pop")

   static std::string DataToHex(const char *data, uint32_t len, bool spaces)
   {
      static constexpr std::array<char, 16> NUM_TO_HEX{{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'}};
      std::string result;
      for (uint32_t i = 0; i < len; i++) {
         result += NUM_TO_HEX[static_cast<uint8_t>(data[i]) >> 4];
         result += NUM_TO_HEX[static_cast<uint8_t>(data[i]) & 0x0f];
         if (spaces && i != len - 1) {
            result += ' ';
         }
      }
      return result;
   }


private:
   // Un-sanitized version of `Load`.
   template <class T>
   static uint64_t __attribute__((no_sanitize_address, no_sanitize_undefined)) DirtyLoad(const T *data)
   {
      static_assert(std::is_same_v<T, uint8_t> || std::is_same_v<T, char>);
      uint64_t result;
      memcpy(&result, data, sizeof(uint64_t));
      return result;
   }
};
// -------------------------------------------------------------------------------------
