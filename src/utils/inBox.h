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

template <typename T>
bool notInBox(const std::array<T, 3> &position, const std::array<T, 3> &low,
              const std::array<T, 3> &high) {
  return not(inBox(position, low, high));
}

} /* namespace autopas */

#endif /* SRC_UTILS_INBOX_H_ */
