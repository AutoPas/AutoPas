/**
 * @file ThreeDimensionalMapping.h
 *
 * @date 15 May 2017
 * @author tchipevn
 */

#pragma once

#include <array>
#include <type_traits>

/**
 * Namespace to handle the conversion between one dimensional and three dimensional indices.
 * The running index is x.
 */
namespace autopas::utils::ThreeDimensionalMapping {

/**
 * Convert a 3d index to a 1d index.
 * @tparam T Type of the indices.
 * @param x x index of the 3d index.
 * @param y y index of the 3d index.
 * @param z z index of the 3d index.
 * @param dims The total dimensions of the index space.
 * @return The 1d index.
 */
template <typename T>
constexpr T threeToOneD(T x, T y, T z, const std::array<T, 3> &dims) {
  static_assert(std::is_integral_v<T>, "threeToOneD requires integral types");
  return (z * dims[1] + y) * dims[0] + x;
}

/**
 * Convert a 3d index to a 1d index.
 * @tparam T Type of the indices.
 * @param index3d The 3d index.
 * @param dims The total dimensions of the index space.
 * @return The 1d index.
 */
template <typename T>
constexpr T threeToOneD(const std::array<T, 3> &index3d, const std::array<T, 3> &dims) {
  static_assert(std::is_integral_v<T>, "threeToOneD requires integral types");
  return (index3d[2] * dims[1] + index3d[1]) * dims[0] + index3d[0];
}

/**
 * Convert a 1d index to a 3d index.
 * @tparam T Type of the indices.
 * @param ind The 1d index.
 * @param dims The total dimensions of the index space.
 * @return The 3d index.
 */
template <typename T>
constexpr std::array<T, 3> oneToThreeD(T ind, const std::array<T, 3> &dims) {
  static_assert(std::is_integral_v<T>, "oneToThreeD requires integral types");
  std::array<T, 3> pos{};
  pos[2] = ind / (dims[0] * dims[1]);
  pos[1] = (ind - pos[2] * dims[0] * dims[1]) / dims[0];
  pos[0] = ind - dims[0] * (pos[1] + dims[1] * pos[2]);
  return pos;
}

// the implementations for the encoding and decoding of Morton indices are taken from this site:
// https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
template<typename T>
constexpr uint64_t threeDtoMortonIndex(T x, T y, T z){
  uint64_t mortonIndex = 0;
  for (uint64_t i = 0; i < (sizeof(uint64_t)* CHAR_BIT)/3; ++i) {
    mortonIndex |= ((x & ((uint64_t)1 << i)) << 2*i) | ((y & ((uint64_t)1 << i)) << (2*i + 1)) | ((z & ((uint64_t)1 << i)) << (2*i + 2));
  }
  return mortonIndex;
}

template<typename T>
constexpr T threeDtoMortonIndex(const std::array<T, 3> &index3d) {
  return threeDtoMortonIndex(index3d[0], index3d[1], index3d[2]);
}

inline uint64_t compactBits(uint64_t n) {
  n &= 0x1249249249249249;
  n = (n ^ (n >> 2)) & 0x30c30c30c30c30c3;
  n = (n ^ (n >> 4)) & 0xf00f00f00f00f00f;
  n = (n ^ (n >> 8)) & 0x00ff0000ff0000ff;
  n = (n ^ (n >> 16)) & 0x00ff00000000ffff;
  n = (n ^ (n >> 32)) & 0x1fffff;
  return n;
}

template<typename T>
constexpr std::array<T, 3> mortonIndexToThreeD(T mortonIndex) {
  std::array<T, 3> pos{};
  pos[2] = compactBits(mortonIndex >> 2);
  pos[1] = compactBits(mortonIndex >> 1);
  pos[0] = compactBits(mortonIndex);
  return pos;
}

}  // namespace autopas::utils::ThreeDimensionalMapping
