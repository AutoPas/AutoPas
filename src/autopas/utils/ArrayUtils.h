/**
 * @file ArrayUtils.h
 *
 * @date 11.07.2019.
 * @author C.Menges
 */

#pragma once

#include <array>
#include <iomanip>
#include <set>
#include <sstream>
#include <vector>

namespace autopas::utils::ArrayUtils {

// specialize a type for containers of type array and vector for use with templated overloaded stream operator.
/**
 * Collection of structs that define what we consider a container.
 */
namespace is_container_impl {
/**
 * Default case: T is not a container.
 * @tparam T
 */
template <typename T>
struct is_container : std::false_type {};
/**
 * Specialization to allow std::array.
 * @tparam T
 * @tparam N
 */
template <typename T, std::size_t N>
struct is_container<std::array<T, N>> : std::true_type {};

/**
 * Specialization to allow std::vector.
 * @tparam Args
 */
template <typename... Args>
struct is_container<std::vector<Args...>> : std::true_type {};

/**
 * Specialization to allow std::vector.
 * @tparam Args
 */
template <typename... Args>
struct is_container<std::set<Args...>> : std::true_type {};

}  // namespace is_container_impl

/**
 * @tparam T Type to check.
 * Type trait to check if a given type  is a container for use with overloaded stream operator.
 * @struct is_container
 * @var is_container::value
 * bool value true if given type is a container false if not
 */
template <typename T>
struct is_container {
  static constexpr bool const value =
      autopas::utils::ArrayUtils::is_container_impl::is_container<std::decay_t<T>>::value;
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
 * Generates a string representation of a container which fulfills the Container requirement (provide cbegin and cend)
 * and appends it to a stream.
 * @tparam Container
 * @param os
 * @param container
 * @param delimiter
 * @param surround
 */
template <class Container>
void to_string(std::ostream &os, const Container &container, const std::string &delimiter = ", ",
               const std::array<std::string, 2> &surround = {"[", "]"}) {
  auto it = std::cbegin(container);
  const auto end = std::cend(container);
  if (it == end) {
    os << surround[0] << surround[1];
    return;
  }
  os << surround[0] << *it;
  for (++it; it != end; ++it) {
    os << delimiter << *it;
  }
  os << surround[1];
}

/**
 * Generates a string representation of a container which fulfills the Container requirement (provide cbegin and cend).
 * @note std::boolalpha is always enabled.
 * @tparam T Type of Container.
 * @param container
 * @param delimiter String that is put between items.
 * @param surround Strings to be put before and after the listing (e.g. brackets).
 * @return String representation of container.
 */
template <class Container>
[[nodiscard]] std::string to_string(const Container &container, const std::string &delimiter = ", ",
                                    const std::array<std::string, 2> &surround = {"[", "]"}) {
  std::ostringstream strStream;
  strStream << std::boolalpha;
  to_string(strStream, container, delimiter, surround);

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

/**
 * Stream operator for containers (array and vector types).
 *
 * This function actually checks if the given Template parameter satisfies is_container.
 * Then Generates a string representation of a container which fulfills the Container requirement (provide cbegin and
 * cend)
 * @tparam array or vector of arbitrary types and sizes
 * @param os string stream
 * @param container
 * @return string representation of a container
 */

template <class Container>
std::enable_if_t<autopas::utils::ArrayUtils::is_container<Container>::value, std::ostream &> operator<<(
    std::ostream &os, const Container &container) {
  const std::string &delimiter = ", ";
  const std::array<std::string, 2> &surround = {"[", "]"};

  to_string(os, container, delimiter, surround);

  return os;
}
}  // namespace autopas::utils::ArrayUtils
