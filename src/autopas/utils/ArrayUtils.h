/**
 * @file ArrayUtils.h
 *
 * @date 11.07.2019.
 * @author C.Menges
 */

#pragma once

#include <array>
#include <sstream>

namespace autopas::utils::ArrayUtils {

/**
 * Creates a new array by performing an element-wise static_cast<>.
 * @tparam output_t Output type.
 * @tparam input_t Input type.
 * @tparam SIZE Size of the array.
 * @param a Input array.
 * @return Array of type std::array<output_t, SIZE>.
 */
template <class output_t, class input_t, std::size_t SIZE>
[[nodiscard]] constexpr std::array<output_t, SIZE> static_cast_array(const std::array<input_t, SIZE> &a) {
  std::array<output_t, SIZE> result{};
  for (std::size_t d = 0; d < SIZE; ++d) {
    result[d] = static_cast<output_t>(a[d]);
  }
  return result;
}

/**
 * Generates a string representation of a container which fulfills the Container requirement (provide cbegin and cend).
 * @note std::boolalpha is always enabled.
 * @tparam T Type of Container.
 * @param a Container.
 * @param delimiter String that is put between items.
 * @param surround Strings to be put before and after the listing (e.g. brackets).
 * @return String representation of a.
 */
template <class T>
[[nodiscard]] std::string to_string(T &&a, const std::string &delimiter = ", ",
                                    const std::array<std::string, 2> &surround = {"[", "]"}) {
  auto it = std::cbegin(a);
  const auto end = std::cend(a);
  if (it == end) {
    return surround[0] + surround[1];
  }
  std::ostringstream strStream;
  strStream << std::boolalpha << surround[0] << *it;
  for (++it; it != end; ++it) {
    strStream << delimiter << *it;
  }
  strStream << surround[1];
  return strStream.str();
}

/**
 * Tests if two vectors (of float-types) are (almost) equal.
 * @tparam T float type
 * @param A vector A
 * @param B vector B
 * @param tol tolerance for equality
 * @return
 */
template <class T>
bool equals(std::vector<T> A, std::vector<T> B, double tol = 1e-10) {
  const auto size = A.size();
  if (size != B.size()) {
    return false;
  }
  for (size_t i = 0; i < size; ++i) {
    if (A[i] - B[i] > tol) {
      return false;
    }
  }
  return true;
}

/**
 * Tests if two vectors of arrays (of float-types) are (almost) equal.
 * @tparam T float type
 * @tparam SIZE size of array
 * @param A vector A
 * @param B vector B
 * @param tol tolerance for equality
 * @return
 */
template <class T, std::size_t SIZE>
bool equals(std::vector<std::array<T, SIZE>> A, std::vector<std::array<T, SIZE>> B, double tol = 1e-10) {
  const auto vec_size = A.size();
  if (vec_size != B.size()) {
    return false;
  }
  for (size_t i = 0; i < vec_size; ++i) {
    for (size_t j = 0; j < SIZE; ++j)
    if (A[i][j] - B[i][j] > tol) {
      return false;
    }
  }
  return true;
}

}  // namespace autopas::utils::ArrayUtils
