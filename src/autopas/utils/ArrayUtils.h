/**
 * @file ArrayUtils.h
 *
 * @date 11.07.2019.
 * @author C.Menges
 */

#pragma once

#include <array>
#include <iomanip>
#include <numeric>
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

/**
 * Given a collection of vectors, redistributes the elements of the vectors so they all have the same (or +1) size.
 * @tparam OuterContainerT Collection type
 * @param vecvec Reference to the collection of vectors to be balanced in place.
 */
template <class OuterContainerT>
void balanceVectors(OuterContainerT &vecvec) {
  balanceVectors(
      vecvec, [](auto &innerContainer) -> auto & { return innerContainer; });
}

template <class OuterContainerT, class F>
void balanceVectors(OuterContainerT &vecvec, F innerContainerToVec) {
  using InnerContainerT = typename OuterContainerT::value_type;
  using ElemT = typename std::remove_reference_t<
      std::invoke_result_t<decltype(innerContainerToVec), InnerContainerT &>>::value_type;
  // calculate vecvec statistics
  const auto vecvecSize = vecvec.size();
  const size_t numElem = std::transform_reduce(vecvec.begin(), vecvec.end(), 0, std::plus<>(),
                                               [&](auto &vec) { return innerContainerToVec(vec).size(); });
  const auto targetSize = static_cast<long>(numElem / vecvecSize);
  auto rest = numElem % vecvecSize;

  std::vector<ElemT> tmpStorage;
  // index of the first subvec that has too few elements
  size_t firstTooFew = 0;
  // repeat as long as there is something in the buffer but at least 2 iterations.
  for (int pass = 0; pass == 0 or not tmpStorage.empty(); ++pass) {
    // scan all subvecs that are not known to have the desired size
    for (size_t i = firstTooFew; i < vecvecSize; ++i) {
      auto &vec = innerContainerToVec(vecvec[i]);
      const auto thisTargetSize = i < rest ? targetSize + 1 : targetSize;
      if (vec.size() > thisTargetSize) {
        // move what is too much to tmpStorage
        const auto startOfTooMuch = vec.begin() + thisTargetSize;
        tmpStorage.insert(tmpStorage.end(), std::make_move_iterator(startOfTooMuch),
                          std::make_move_iterator(vec.end()));
        vec.erase(startOfTooMuch, vec.end());
        // if firstTooFew points here bump it, since this vector is now satisfied
        if (firstTooFew == i) {
          ++firstTooFew;
        }
      } else if (vec.size() < thisTargetSize) {
        // fill too small vectors from storage
        vec.reserve(thisTargetSize);
        const long numMissing = thisTargetSize - vec.size();
        const auto startOfInsertion =
            tmpStorage.begin() + (std::max(0l, static_cast<long>(tmpStorage.size()) - numMissing));
        vec.insert(vec.end(), std::make_move_iterator(startOfInsertion), std::make_move_iterator(tmpStorage.end()));
        tmpStorage.erase(startOfInsertion, tmpStorage.end());

      } else if (firstTooFew == i) {
        // if firstTooFew points here bump it, since this vector is already satisfied
        ++firstTooFew;
      }
      // after the first pass all excess elements are collected. So as soon as everything is distributed we are done.
      if (pass > 0 and tmpStorage.empty()) {
        break;
      }
    }
  }
}

}  // namespace autopas::utils::ArrayUtils
