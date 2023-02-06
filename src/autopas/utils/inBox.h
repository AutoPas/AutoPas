/**
 * @file inBox.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <type_traits>

namespace autopas::utils {
/**
 * Checks if position is inside of a box defined by low and high.
 * The lower corner is included in, the upper is excluded from the box.
 * i.e. [low[0], high[0]) x [low[1], high[1]) x [low[2], high[2])
 *
 * @tparam T the type of floating point check
 * @param position the position that should be checked
 * @param low the lower corner of the box (inclusive)
 * @param high the upper corner of the box (exclusive)
 * @return true if position is inside the box, false otherwise
 */
template <typename T>
bool inBox(const std::array<T, 3> &position, const std::array<T, 3> &low, const std::array<T, 3> &high) {
  static_assert(std::is_floating_point<T>::value, "inBox assumes floating point types");

  bool inBox = true;
  for (int d = 0; d < 3; ++d) {
    const bool isLargerThanLower = position[d] >= low[d];
    const bool isSmallerThanHigher = position[d] < high[d];
    inBox = inBox and isLargerThanLower and isSmallerThanHigher;
  }
  return inBox;
}

/**
 * Checks if position is not inside of a box defined by low and high.
 * The lower corner is included in, the upper is excluded from the box.
 * i.e. [low[0], high[0]) x [low[1], high[1]) x [low[2], high[2])
 *
 * @tparam T the type of floating point check
 * @param position the position that should be checked
 * @param low the lower corner of the box (inclusive)
 * @param high the upper corner of the box (exclusive)
 * @return true if position is not inside the box, false otherwise
 */
template <typename T>
bool notInBox(const std::array<T, 3> &position, const std::array<T, 3> &low, const std::array<T, 3> &high) {
  return not(inBox(position, low, high));
}

/**
 * Checks if two boxes have overlap.
 * @tparam T
 * @param boxALow
 * @param boxAHigh
 * @param boxBLow
 * @param boxBHigh
 * @return
 */
template <typename T>
bool boxesOverlap(const std::array<T, 3> &boxALow, const std::array<T, 3> &boxAHigh, const std::array<T, 3> &boxBLow,
                  const std::array<T, 3> &boxBHigh) {
  static_assert(std::is_floating_point<T>::value, "boxesOverlap assumes floating point types");

  auto overlap1D = [&](size_t dim) { return boxAHigh[dim] >= boxBLow[dim] and boxBHigh[dim] >= boxALow[dim]; };

  return overlap1D(0) and overlap1D(1) and overlap1D(2);
}
}  // namespace autopas::utils
