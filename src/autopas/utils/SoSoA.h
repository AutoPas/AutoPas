/**
 * @file SoSoA.h
 * @authors S. J. Newcome
 * @date 03/09/2024
 */

#pragma once

#include <algorithm>
#include <initializer_list>
#include <map>
#include <tuple>
#include <vector>

#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/SoAStorage.h"
#include "autopas/utils/SoAType.h"
#include "autopas/utils/SoSoAView.h"

namespace autopas {

/**
 * Structure of Structure of Arrays class. This is a tuple of vectors of SoAs of different types.
 *
 * It is intended for use with multi-component particles: where each particle contributes to different types of SoAs
 * representing different types of components and each particle can contribute to zero, one, or more SoAs of a certain
 * type depending on the number of components of that type.
 *
 * For each SoAType: there exists a vector of SoAs. We call the number of SoAs of a certain type the maximum "depth" of
 * that type. The depth must be larger than or equal to the maximum number of components of that type of any molecule in the
 * SoSoA. We consider the SoSoA to have a "length" equal to the number entries in each array i.e. the number of particles. This
 * should always be equal for every array.
 *
 * If a particle does not have as many components of a type as the depth of that type, "null" entries are used.
 *
 * For example, an SoSoA containing 20 TIP4 multi-site water molecules consisting of 1 Lennard-Jones site and 3
 * electrostatic sites. A (well-fitting) SoSoA would contain a maximum depth of 1 Lennard-Jones-Type SoAs and a maximum depth of 3
 * Electrostatic-Type SoAs, with all arrays of length 20. If the SoSoA also contained a molecule consisting of 2 Lennard
 * -Jones sites and 0 electrostatic site, the SoSoA would contain a maximum depth of 2 LJ-Type SoAs and a maximum depth of 3 Electrostatic-Type
 * SoAs, with null entries in the deeper SoAs.
 *
 * @warning This way of handling SoAs for multi-component particles is intended primarily for systems of only a few
 * different particle types. It is probably fairly inefficient for multiple particle types with different numbers of
 * components -> but could be made less so with sorting such that particles are ordered in descending number of
 * components of each type.
 *
 * @tparam SoAtypes List of  SoAType structs that make up this SoSoA
 */
template <typename... SoAtypes>
class SoSoA {
 public:
  /**
   * Default constructor.
   */
  SoSoA() = default;

  /**
   * Copy constructor.
   * @param sosoa SoSoA to copy.
   */
  SoSoA(const SoSoA &sosoa) = default;

  using SoSoAType = std::tuple<std::vector<SoA<SoAtypes>...>>;

  /**
   * Number of types
   */
  static const size_t numTypes = std::tuple_size<SoSoAType>::value;

  /**
   * The actual SoSoA.
   */
  SoSoAType SoAs;

  /**
   * Resizes all SoAs to the given length.
   * @param length new length.
   */
  void resizeLengthOfArrays(size_t length) {
    resizeLengthArraysImpl(length, std::make_index_sequence<numTypes>{});
  }

  /**
   * Appends the other SoSoA buffer to this in a lengthwise manner.
   * @param other other SoSoA buffer
   */
  void append(const SoSoAType &other) {
    if (other.size() > 0) {
      appendImpl(other, std::make_index_sequence<numTypes>{});
    }
  }

  /**
   * Appends the other SoSoA buffer to this in a lengthwise manner.
   * @param other other SoSoA buffer as an SoSoAView
   */
  void append(const SoSoAView<SoSoAType> &other) {
    if (other.size() > 0) {
      appendImpl(other, std::make_index_sequence<numTypes>{});
    }
  }

  /**
   * Returns a pointer to the given attribute vector in the given SoA.
   * @tparam attribute
   * @param soaTypeIndex type index of given SoA
   * @param depth
   * @return
   */
  template <size_t attribute>
  auto begin(size_t soaTypeIndex, size_t depth) {
    return SoAs[soaTypeIndex][depth].template begin<attribute>();
  }

  /**
   * Returns the "length" of the SoSoA i.e. the number of particles. This is the maximum length of any SoA because a
   * SoA could have less entries than the number of particles.
   *
   * @return Number of particles.
   */
  inline size_t size() const {
    size_t maxLength = 0;
    utils::TupleUtils::for_each(SoAs.getTuple(), [&](auto &SoAsOfSingleType) {
      for (auto &soa : SoAsOfSingleType) {
        maxLength = std::max(maxLength, soa->size());
      }
    });
    return maxLength;
  }

  /**
   * Delete all particles in the SoSoA
   */
  void clear() {
    utils::TupleUtils::for_each(SoAs.getTuple(), [&](auto &SoAsOfSingleType) {
      for (auto &soa : SoAsOfSingleType) {
        soa->clear();
      }
    });
  }

  /**
   * Constructs a SoSoAView for the whole SoSoA and returns it.
   * @return the constructed SoSoAView on the whole SoSoA.
   */
  SoSoAView<SoSoAType> constructView() { return {this, 0, size()}; }

  /**
   * Constructs a view that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `this->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param startIndex The index of the first entry to view.
   * @param endIndex The index of the entry after the last entry to view.
   * @return the constructed SoSoAView from \p startIndex (inclusive) to \p endIndex (exclusive).
   */
  SoSoAView<SoSoAType> constructView(size_t startIndex, size_t endIndex) { return {this, startIndex, endIndex}; }

 private:
  /**
   * Helper for clarity to get the maximum "depth" of a given SoA type i.e. how many SoAs are there of this type.
   * @param soaTypeIndex index of the type of SoA.
   * @return maximum "depth" of the given SoA type.
   */
  size_t getMaxDepthOfGivenType(size_t soaTypeIndex) {
    return SoAs[soaTypeIndex].size();
  }

  /**
   * Helper for clarity to resize the maximum "depth" of a given SoA type.
   * @param soaTypeIndex index of the type of SoA.
   * @param newMaxDepth new maximum depth the vector of SoAs is resized to.
   */
  void resizeMaxDepthOfGivenType(size_t soaTypeIndex, size_t newMaxDepth) {
    SoAs[soaTypeIndex].resize(newMaxDepth);
  }

  /**
   * Helper for the actual implementation of resizeLengthOfArrays.
   * @tparam soaTypeIndices sequence of size_t indices representing the SoA types.
   * @param length length arrays are resized to.
   */
  template <size_t... soaTypeIndices>
  void resizeLengthOfArraysImpl(size_t length, std::index_sequence<soaTypeIndices...>) {
    // fold expression
    (resizeLengthOfArraysOfSingleType<soaTypeIndices>(length), ...);
  }

  /**
   * Applies resizeArrays to every SoA of a given type.
   * @tparam soaTypeIndex soaTypeIndex type index of given SoA
   * @param length length arrays are resized to.
   */
  template <size_t soaTypeIndex>
  void resizeLengthOfArraysOfSingleType(size_t length) {
    for (auto soa : SoAs[soaTypeIndex] ) {
      soa.resizeArrays(length);
    }
  }

  /**
   * Helper for the actual implementation of append.
   * @tparam soaTypeIndices sequence of size_t indices representing the SoA types.
   * @param other other SoSoA buffer.
   */
  template <size_t... soaTypeIndices>
  void appendImpl(SoSoAType &other, std::index_sequence<soaTypeIndices...>) {
    // fold expression
    (appendSingleType<soaTypeIndices>(other), ...);
  }

  /**
   * Helper for the actual implementation of append.
   * @tparam soaTypeIndices sequence of size_t indices representing the SoA types.
   * @param other other SoSoAView buffer.
   */
  template <size_t... soaTypeIndices>
  void appendImpl(SoSoAView<SoSoAType> &other, std::index_sequence<soaTypeIndices...>) {
    // fold expression
    (appendSingleType<soaTypeIndices>(other), ...);
  }

  /**
   * Applies append to every SoA of a given type with the relevant SoA from the other SoSoA buffer. Increases the maximum
   * depth of the SoA type if the other SoSoA buffer has deeper SoAs of this type.
   * @tparam soaTypeIndex soaTypeIndex type index of given SoA
   * @param other other SoSoA buffer.
   */
  template <size_t soaTypeIndex>
  void appendSingleType(SoSoAType &other) {
    // check if resize is needed and if so resize
    const auto otherDepth = other[soaTypeIndex].size();
    if (otherDepth > SoAs[soaTypeIndex].size()) {
      SoAs[soaTypeIndex].resize(otherDepth);
    }

    // Append SoAs up until otherDepth.
    for (size_t i = 0; i < otherDepth; ++i) {
      SoAs[soaTypeIndex][i].append(other[soaTypeIndex][i]);
    }
  }

  /**
   * Applies append to every SoA of a given type with the relevant SoAView from the other SoSoAView buffer. Increases the maximum
   * depth of the SoA type if the other SoSoAView buffer has deeper SoAs of this type.
   * @tparam soaTypeIndex soaTypeIndex type index of given SoA
   * @param other other SoSoAView buffer.
   */
  template <size_t soaTypeIndex>
  void appendSingleType(SoSoAView<SoSoAType> &other) {
    // check if resize is needed and if so resize
    const auto otherDepth = other[soaTypeIndex].size();
    if (otherDepth > SoAs[soaTypeIndex].size()) {
      SoAs[soaTypeIndex].resize(otherDepth);
    }

    // Append SoAs up until otherDepth.
    for (size_t i = 0; i < otherDepth; ++i) {
      // A SoSoAView is a view on a SoSoA, not a Structure of SoAViews. Therefore, the SoAViews must be constructed from
      // the SoAs.
      auto soaView = SoAView(other[soaTypeIndex][i], other->_s);
      SoAs[soaTypeIndex][i].append(soaView);
    }
  }
};
}  // namespace autopas
