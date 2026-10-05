/**
 * @file HashCombine.h
 * @date 14/09/2026
 * @author S. J. Newcome
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>

namespace autopas::utils {

/**
 * Hash of a single value.
 *
 * Uses std::hash if it is specialized for T (arithmetic types, enums, strings, pointers, ...) and otherwise falls back
 * to a static_cast to std::size_t. The latter covers AutoPas Option types, which implicitly convert to their underlying
 * enum.
 *
 * @tparam T Type of the value.
 * @param value Value to hash.
 * @return Hash of value.
 */
template <class T>
std::size_t hashValue(const T &value) {
  if constexpr (requires { std::hash<T>{}(value); }) {
    return std::hash<T>{}(value);
  } else {
    return static_cast<std::size_t>(value);
  }
}

/**
 * Combines the hashes of an arbitrary number of values into one hash.
 *
 * Equivalent to repeated application of boost::hash_combine (Boost >= 1.81), i.e.
 * `seed = mix(seed + 0x9e3779b9 + hashValue(value))`.
 *
 * The mixer is the xmxmx bit mixer used by boost::hash_combine (boost/container_hash/detail/hash_mix.hpp). It is
 * reproduced here instead of pulling in Boost, and as we use this for ease and not in performance critical situations,
 * we do not care so much about getting the latest version from Boost.
 *
 * The parameters used are those from Jon Maiga's implementation (https://jonkagstrom.com/mx3/mx3_rev2.html).
 *
 * @tparam Values_T Types of the values to hash.
 * @param values Values to hash.
 * @return Combined hash.
 */
template <class... Values_T>
std::size_t hashCombine(const Values_T &...values) {
  static_assert(sizeof(std::size_t) == 8,
                "hashCombine() only supports a 64-bit std::size_t, as the mixing constants are chosen for that width.");

  constexpr auto mix = [](std::size_t x) {
    constexpr std::uint64_t m = 0xe9846af9b1a615dULL;
    x ^= x >> 32;
    x *= m;
    x ^= x >> 32;
    x *= m;
    x ^= x >> 28;
    return x;
  };

  std::size_t seed = 0;
  ((seed = mix(seed + 0x9e3779b9 + hashValue(values))), ...);
  return seed;
}

}  // namespace autopas::utils