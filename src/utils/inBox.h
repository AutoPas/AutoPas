/*
 * inBox.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_UTILS_INBOX_H_
#define SRC_UTILS_INBOX_H_

#include <array>
#include <type_traits>

namespace autopas {

/**
 * Checks if position is inside of a box defined by low and high.
 * The lower corner is include, the upper is exclusive to the box.
 *
 * @tparam T the type of floating point check
 * @param position the position that should be checked
 * @param low the lower corner of the box
 * @param high the upper corner of the box
 * @return true if position is inside the box, false otherwise
 */
template <typename T>
bool inBox(const std::array<T, 3> &position, const std::array<T, 3> &low,
           const std::array<T, 3> &high) {
  static_assert(std::is_floating_point<T>::value,
                "inBox assumes floating point types");

  bool inBox = true;
  for (int d = 0; d < 3; ++d) {
    const bool isLargerThanLower = position[d] >= low[d];
    const bool isSmallerThanHigher = position[d] < high[d];
    inBox &= isLargerThanLower and isSmallerThanHigher;
  }
  return inBox;
}

/**
 * Checks if position is not inside of a box defined by low and high.
 * The lower corner is include, the upper is exclusive to the box.
 *
 * @tparam T the type of floating point check
 * @param position the position that should be checked
 * @param low the lower corner of the box
 * @param high the upper corner of the box
 * @return true if position is not inside the box, false otherwise
 */
template <typename T>
bool notInBox(const std::array<T, 3> &position, const std::array<T, 3> &low,
              const std::array<T, 3> &high) {
  return not(inBox(position, low, high));
}

} /* namespace autopas */

#endif /* SRC_UTILS_INBOX_H_ */
