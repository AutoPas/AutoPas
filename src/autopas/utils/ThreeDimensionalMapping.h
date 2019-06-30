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
  static_assert(std::is_integral<T>::value, "threeToOneD requires integral types");
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
  static_assert(std::is_integral<T>::value, "threeToOneD requires integral types");
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
  static_assert(std::is_integral<T>::value, "oneToThreeD requires integral types");
  std::array<T, 3> pos{};
  pos[2] = ind / (dims[0] * dims[1]);
  pos[1] = (ind - pos[2] * dims[0] * dims[1]) / dims[0];
  pos[0] = ind - dims[0] * (pos[1] + dims[1] * pos[2]);
  return pos;
}

}  // namespace autopas::utils::ThreeDimensionalMapping
