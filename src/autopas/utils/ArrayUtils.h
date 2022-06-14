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

// specialize a type for all of the STL containers.
/**
 * Collection of structs that define what we consider a container. Remove / add
 * whatever you need.
 */
namespace is_container_impl {
/**
 * Default case: T is not a container.
 * @tparam T
 */
template <typename T> struct is_container : std::false_type {};
/**
 * Specialization to allow std::array.
 * @tparam T
 * @tparam N
 */
template <typename T, std::size_t N>
struct is_container<std::array<T, N>> : std::true_type {};

} // namespace is_container_impl

/**
 * Type trait to check if a given type is a container.
 * @tparam T Type to check.
 */
template <typename T> struct is_container {
  static constexpr bool const value =
      is_container_impl::is_container<std::decay_t<T>>::value;
};

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
    strStream << delimiter << *it;
  }
  strStream << surround[1];
  return strStream.str();
}

}  // namespace autopas::utils::ArrayUtils

/**
 * Stream operator for containers.
 *
 * This function actually checks if the given Template parameter satisfies is_container.
 *
 * @tparam Container
 * @param os
 * @param container
 * @return
 */
template <class Container>
std::enable_if_t<autopas::utils::ArrayUtils::is_container<Container>::value, std::ostream &>
operator<<(std::ostream &os, const Container &container) {
  os << autopas::utils::ArrayUtils::to_string(container);
  return os;
}
