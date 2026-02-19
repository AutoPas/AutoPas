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
#include <sstream>
#include <vector>

#include "autopas/utils/ContainerConcept.h"

namespace autopas::utils::ArrayUtils {
/**
 * Creates a new array by performing an element-wise static_cast<>.
 *
 * @note This function returns a new copy of the array with the desired type!
 *
 * Even though this is implemented to copy the array, compilers optimize this away: https://gcc.godbolt.org/z/6dav1PEGP
 *
 * @tparam output_t Output type.
 * @tparam input_t Input type.
 * @tparam SIZE Size of the array.
 * @param a Input array.
 * @return Array of type std::array<output_t, SIZE>.
 */
template <class output_t, class input_t, std::size_t SIZE>
[[nodiscard]] constexpr std::array<output_t, SIZE> static_cast_copy_array(const std::array<input_t, SIZE> &a) {
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
 * @tparam Fun Function type (Container::element) -> implicit std::string
 * @param os
 * @param container
 * @param elemToString Function converting one element of container to something that is implicitly convertible to
 * std::string.
 * @param delimiter
 * @param surround
 */
template <class Container, class Fun>
void to_string(std::ostream &os, const Container &container, const std::string &delimiter,
               const std::array<std::string, 2> &surround, Fun elemToString) {
  auto it = std::cbegin(container);
  const auto end = std::cend(container);
  if (it == end) {
    os << surround[0] << surround[1];
    return;
  }
  os << surround[0] << elemToString(*it);
  for (++it; it != end; ++it) {
    os << delimiter << elemToString(*it);
  }
  os << surround[1];
}

/**
 * Version of to_string() with simpler signature and default arguments.
 * @tparam Container
 * @param os
 * @param container
 * @param delimiter
 * @param surround
 */
template <class Container>
void to_string(std::ostream &os, const Container &container, const std::string &delimiter = ", ",
               const std::array<std::string, 2> &surround = {"[", "]"}) {
  to_string(os, container, delimiter, surround, [](const auto &x) { return x; });
}

/**
 * Generates a string representation of a container which fulfills the Container requirement (provide cbegin and cend).
 * @note std::boolalpha is always enabled.
 * @tparam T Type of Container.
 * @tparam Fun Function type (Container::element) -> implicit std::string
 * @param container
 * @param elemToString Function converting one element of container to something that is implicitly convertible to
 * std::string.
 * @param delimiter String that is put between items.
 * @param surround Strings to be put before and after the listing (e.g. brackets).
 * @return String representation of container.
 */
template <class Container, class Fun>
[[nodiscard]] std::string to_string(const Container &container, const std::string &delimiter,
                                    const std::array<std::string, 2> &surround, Fun elemToString) {
  std::ostringstream strStream;
  strStream << std::boolalpha;
  to_string(strStream, container, delimiter, surround, elemToString);

  return strStream.str();
}

/**
 * Version of to_string() with simpler signature and default arguments.
 * @tparam Container
 * @param container
 * @param delimiter
 * @param surround
 * @return
 */
template <class Container>
[[nodiscard]] std::string to_string(const Container &container, const std::string &delimiter = ", ",
                                    const std::array<std::string, 2> &surround = {"[", "]"}) {
  return to_string(container, delimiter, surround, [](const auto &elem) { return elem; });
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

template <ContainerType Container>
std::ostream &operator<<(std::ostream &os, const Container &container) {
  const std::string &delimiter = ", ";
  const std::array<std::string, 2> &surround = {"[", "]"};

  to_string(os, container, delimiter, surround);

  return os;
}

/**
 * Given a collection of vectors, redistributes the elements of the vectors so they all have the same (or +1) size.
 *
 * Elements are taken from the ends of too-long vectors and appended to the ends of too-short vectors.
 * The overall ordering of elements is not preserved.
 *
 * @tparam OuterContainerT Collection type
 * @param vecvec Reference to the collection of vectors to be balanced in place.
 */
template <class OuterContainerT>
void balanceVectors(OuterContainerT &vecvec) {
  balanceVectors(
      vecvec, [](auto &innerContainer) -> auto & { return innerContainer; });
}

/**
 * Given a collection of containers that hold vectors,
 * redistributes the elements of the vectors so they all have the same (or +1) size.
 *
 * Elements are taken from the ends of too-long vectors and appended to the ends of too-short vectors.
 * The overall ordering of elements is not preserved.
 *
 * @tparam OuterContainerT Collection type
 * @tparam F Type of the function innerContainerToVec
 * @param vecvec Reference to the collection of vectors to be balanced in place.
 * @param innerContainerToVec Function to map inner containers to std::vector&.
 */
template <class OuterContainerT, class F>
void balanceVectors(OuterContainerT &vecvec, F innerContainerToVec) {
  using InnerContainerT = typename OuterContainerT::value_type;
  using ElemT = typename std::remove_reference_t<
      std::invoke_result_t<decltype(innerContainerToVec), InnerContainerT &>>::value_type;
  // calculate vecvec statistics
  const auto vecvecSize = vecvec.size();
  if (vecvecSize == 0) {
    // nothing to do
    return;
  }
  const size_t numElem = std::transform_reduce(vecvec.begin(), vecvec.end(), 0, std::plus<>(),
                                               [&](auto &vec) { return innerContainerToVec(vec).size(); });
  const auto targetSize = static_cast<long>(numElem / vecvecSize);
  const auto remainder = numElem % vecvecSize;

  std::vector<ElemT> tmpStorage;
  // index of the first subvec that has too few elements
  size_t firstTooFew = 0;
  // repeat as long as there is something in the buffer but at least 2 iterations.
  for (int pass = 0; pass == 0 or not tmpStorage.empty(); ++pass) {
    // scan all subvecs that are not known to have the desired size
    for (size_t i = firstTooFew; i < vecvecSize; ++i) {
      auto &vec = innerContainerToVec(vecvec[i]);
      const auto thisTargetSize = i < remainder ? targetSize + 1 : targetSize;
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
