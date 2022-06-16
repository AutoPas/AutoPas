/**
 * @file ArrayUtils.h
 *
 * @date 11.07.2019.
 * @author C.Menges
 */

#pragma once

#include <array>
#include <sstream>
#include <iomanip>  

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
 * @return String representation of container.
 */
template <class Container>
[[nodiscard]] std::string
to_string(const Container &container, const std::string &delimiter = ", ", const std::array<std::string, 2> &surround = {"[", "]"}) {
  auto it = std::cbegin(container);
  const auto end = std::cend(container);
  if (it == end) {
    return surround[0] + surround[1];
  }
  std::ostringstream strStream;
  strStream << std::boolalpha << surround[0] << *it;
  for (++it; it != end; ++it) {
    strStream << delimiter << std::setprecision(2)<<*it;
  }
  strStream << surround[1];
  return strStream.str();
}

}  // namespace autopas::utils::ArrayUtils

