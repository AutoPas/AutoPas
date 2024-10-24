/**
 * @file SoA.h
 * @authors S. J. Newcome
 * @date 03/09/2024
 */

#pragma once

#include <tuple>
#include <vector>

#include "autopas/utils/SoAPartition.h"
#include "autopas/utils/SoAPartitionView.h"
#include "autopas/utils/SoAView.h"
#include "autopas/utils/TupleUtils.h"

namespace autopas {

/**
 * Structure of Arrays class. This is a pair of a main SoA partition and a tuple of zero, one, or more different types of vectors
 * of SoA partitions (See SoAType.h)
 *
 * For any particle with a fixed number of attributes, only the main SoA partition is to be used. This includes
 * particle attributes of a fixed number per particle (e.g. there are always 3 position buffers per particle).
 *
 * The other partitions are intended for use with multi-component particles: where each particle may have a variable
 * number of components (not known at compile time) and these components may come in a variety of types (known at
 * compile type).
 *
 * This is implemented as follows:
 * - The SoA has a "length" equal to the number of particles.
 * - The main SoA partition consists of an array of that length for each attribute of that particle (attributes are
 *   defined in the particle class)
 * - For each additional component type, there exists a vector of SoA partitions. We call the number of partitions (i.e.
 *   the size of the vector), the maximum "depth" of that type of SoA partition. This should be larger than or equal to
 *   the maximum number of components of that type of any particle in the SoA.
 * - If a particle does not have as many components of a type as the depth of that type, "null" entries are used.
 *
 * Illustrative Example
 *
 * Consider a SoA containing 20 TIP4 multi-site water molecules consisting of 1 Lennard-Jones site and 3 electrostatic
 * sites. A (well-fitting) SoA would contain
 * - a main partition consisting of molecule-level attributes (e.g. position of center of mass, total force, ownership)
 * - 1 Lennard-Jones-Type partitions containing LJ-attributes (e.g. site position, epsilon, sigma) (i.e. max depth 1)
 * - 3 Electrostatic-Type partitions (e.g. containing site position, charge) (i.e. max depth 3)
 * with all arrays of length 20. The molecular model could also allow for a third type of site (e.g. dipole) that is not
 * used i.e. there exists a depth of 0 and no memory is allocated for this.
 *
 * If the SoA also contained a molecule consisting of 2 Lennard-Jones sites and 0 electrostatic site, the SoA would
 * contain a maximum depth of 2 LJ-Type partitions and a maximum depth of 3 Electrostatic-Type partitions, with null
 * entries in the deeper partitions.
 *
 * @warning This way of handling SoAs for multi-component particles is intended primarily for systems of only a few
 * different particle types. It is probably fairly inefficient for multiple particle types with different numbers of
 * components -> but could be made less so with sorting such that particles are ordered in descending number of
 * components of each type, and deeper partitions are resized to tightly fit the number of sites at that depth.
 *
 * @tparam MainPartitionType The type of the main SoA partition
 * @tparam AdditionalPartitionTypes List of addition partition types that make up this SoA
 */
template <class SoAType>
class SoA {
 public:
  /**
   * Default constructor.
   */
  SoA() = default;

  /**
   * Copy constructor.
   * @param SoA SoA to copy.
   */
  SoA(const SoA &SoA) = default;

  /**
   * Type alias for SoAView of this type.
   */
  using SoAViewType = SoAView<SoAType>;

  /**
   * The Main SoA Partition Type
   */
  using MainPartitionType = typename SoAType::MainPartitionType;

  /**
   * The Additional SoA Partitions Type. This is a std::tuple<std::vector<AdditionalSoAPartitionTypes>...> (a tuple,
   * with an element per partition type, of vectors of partitions of that type (which themselves are tuples of vectors)).
   */
  using AdditionalPartitionType = typename SoAType::AdditionalPartitionsType;

  using AdditionalPartitionTypesReference = typename SoAType::AdditionalPartitionTypesReference;

  /**
   * Number of additional partition types
   */
  static const size_t numAdditionalTypes = std::tuple_size<AdditionalPartitionType>::value;

  /**
   * Resizes all partitions to the given length.
   * @param length new length.
   */
  void resizeLengthOfArrays(size_t length) {
    _mainSoAPartition.resizeArrays(length);
    resizeLengthOfArraysAdditionalTypesImpl(length, std::make_index_sequence<numAdditionalTypes>{});
  }

  /**
   * Appends the other SoA buffer to this in a lengthwise manner.
   * @param other other SoA buffer
   */
  void append(const SoA<SoAType> &other) {
    _mainSoAPartition.append(other._mainSoAPartition);
    if (other.size() > 0) {
      appendAdditionalTypesImpl(other, std::make_index_sequence<numAdditionalTypes>{});
    }
  }

  /**
   * Appends the other SoA buffer to this in a lengthwise manner.
   * @param other other SoA buffer as an SoAView
   */
  void append(const SoAViewType &other) {
    _mainSoAPartition.append(other._mainSoAPartitionView);
    if (other.size() > 0) {
      appendAdditionalTypesImpl(other, std::make_index_sequence<numAdditionalTypes>{});
    }
  }
  
  /**
   * Returns a pointer to the given attribute vector in the main SoA partition.
   * @tparam attribute 
   * @return 
   */
  template <size_t attribute>
  auto begin() {
    return _mainSoAPartition.template begin<attribute>();
  }

  /**
   * Returns a pointer to the given attribute vector in the given additional SoA partition (of the given type and depth).
   * @tparam additionalPartitionTypeIndex index corresponding to the desired additional SoA partition's type
   * @tparam attribute desired attribute index
   * @param depth depth of desired SoA partition
   * @return
   */
  template <size_t additionalPartitionTypeIndex, size_t attribute>
  auto begin(size_t depth) {
    return std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions)[depth].template begin<attribute>();
  }

  /**
   * Returns pointers to the start of all vectors of the given desired attributes in both main and additional partitions.
   *
   * For the main partition, pointers are simply a tuple of pointers to the beginning of arrays corresponding to
   * the given mainAttributes.
   *
   * For the additional partitions, pointers are arranged in a structure that mimics the SoA's additional partitions
   * structure:
   *                                       (A)      (B)        (C)       (D)
   * SoA additional partition structure: tuple of vectors of tuples of vectors
   * Returned pointers:                  tuple of vectors of tuples of pointers (to start of vectors D)
   *
   * Note: For returned pointers, the tuple C is the requested subset of the tuple C of the structure. The sizes of the
   * vectors B are equal.
   *
   * @tparam mainAttributesType the type of the array of desired mainAttributes.
   * @tparam mainAttributes the array of desired mainAttributes.
   * @tparam additionalAttributesType the type of the tuple of arrays of desired additionalAttributes. There should be
   * one tuple element for each type of additional SoA partition type and, for each type, one array element for each desired
   * attribute of that type.
   * @tparam additionalAttr the desired additionalAttributes.
   * @return a pair of
   * * a tuple (an element for each desired mainAttribute) of pointers to the start of the corresponding arrays
   * * a tuple (an element for each SoAPartitionType) of a vector (of size maxDepth for each SoAPartitionType) of
   *   a tuple (an element for each desired additionalAttributeof that SoAPartitionType) of pointers to the start of the
   *   corresponding arrays.
   */
  template <typename Functor>
  auto begin() {
    static constexpr auto numMainAttr = Functor::getNeededAttr().size();
    return std::make_pair(beginMainAttrImpl<Functor>(std::make_index_sequence<numMainAttr>{}),
        beginAdditionalAttrImpl<Functor>(std::make_index_sequence<numAdditionalTypes>{})
    );

  }

  /**
   * Returns the "length" of the SoA i.e. the number of particles. This is the size of the main partition.
   *
   * @return Number of particles.
   */
  [[nodiscard]] inline size_t size() const {
    return _mainSoAPartition.size();
  }

  /**
   * Delete all particles in the SoA
   */
  void clear() {
    _mainSoAPartition.clear();
    utils::TupleUtils::for_each(_additionalSoAPartitions.getTuple(), [&](auto &SoAsOfSingleType) {
      for (auto &soaTmp : SoAsOfSingleType) {
        soaTmp->clear();
      }
    });
  }

  /**
   * Constructs a SoAView for the whole SoA and returns it.
   * @return the constructed SoAView on the whole SoA.
   */
  SoAViewType constructView() { return {this, 0, size()}; }

  /**
   * Constructs a view that starts at \p startIndex (inclusive) and ends at \p endIndex (exclusive).
   *
   * \p startIndex and \p endIndex have to be between 0 (inclusive) and `this->size()` (inclusive). \p
   * endIndex has to be greater or equal to \p startIndex.
   * @param startIndex The index of the first entry to view.
   * @param endIndex The index of the entry after the last entry to view.
   * @return the constructed SoAView from \p startIndex (inclusive) to \p endIndex (exclusive).
   */
  SoAViewType constructView(size_t startIndex, size_t endIndex) {
    return {this, startIndex, endIndex}; 
  }

  /**
   * Helper for clarity to get the maximum "depth" of a given additional partition type i.e. how many partitions are
   * there of this type.
   * @tparam additionalPartitionTypeIndex index corresponding to the desired additional SoA partition's type.
   * @return maximum "depth" of the given partition type.
   */
  template <size_t additionalPartitionTypeIndex>
  size_t getMaxDepthOfGivenType() {
    return std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions).size();
  }

  /**
   * Helper for clarity to get the maximum "depths" of all additional partition types i.e. how many partitions are there
   * of each type.
   * @return array of max depths
   */
  std::array<size_t, numAdditionalTypes> getMaxDepths() {
    return getMaxDepthsImpl(std::make_index_sequence<numAdditionalTypes>{});
  }

  /**
   * Helper for clarity to resize the maximum "depth" of a given additional partition type.
   * @tparam additionalPartitionTypeIndex index corresponding to the desired additional SoA partition's type.
   * @param newMaxDepth new maximum depth the vector of partitions is resized to.
   */
  template <size_t additionalPartitionTypeIndex>
  void resizeMaxDepthOfGivenType(size_t newMaxDepth) {
    std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions).resize(newMaxDepth);
  }

  /**
   * The actual SoA data:
   */

  /**
   * The Main SoA Partition.
   */
  MainPartitionType _mainSoAPartition;

  /**
   * The Additional SoA Partitions.
   */
  AdditionalPartitionType _additionalSoAPartitions;

 private:
  /**
   * Helper for the actual implementation of resizeLengthOfArrays that handles additional partition types.
   * @tparam additionalPartitionTypeIndices sequence of size_t indices representing the partition types.
   * @param length length arrays are resized to.
   */
  template <size_t... additionalPartitionTypeIndices>
  void resizeLengthOfArraysAdditionalTypesImpl(size_t length, std::index_sequence<additionalPartitionTypeIndices...>) {
    // fold expression
    (resizeLengthOfArraysOfSingleType<additionalPartitionTypeIndices>(length), ...);
  }

  /**
   * Applies resizeArrays to every additional partition of a given type.
   * @tparam additionalPartitionTypeIndex index corresponding to the desired additional SoA partition's type.
   * @param length length arrays are resized to.
   */
  template <size_t additionalPartitionTypeIndex>
  void resizeLengthOfArraysOfSingleType(size_t length) {
    for (auto soaTmp : std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions) ) {
      soaTmp.resizeArrays(length);
    }
  }

  /**
   * Helper for the actual implementation of append that handles additional partition types.
   * @tparam additionalPartitionTypeIndices sequence of size_t indices representing the partition types.
   * @param other other SoA buffer.
   */
  template <size_t... additionalPartitionTypeIndices>
  void appendAdditionalTypesImpl(const SoA<SoAType> &other, std::index_sequence<additionalPartitionTypeIndices...>) {
    // fold expression
    (appendSingleType<additionalPartitionTypeIndices>(other), ...);
  }

  /**
   * Helper for the actual implementation of append that handles additional partition types.
   * @tparam additionalPartitionTypeIndices sequence of size_t indices representing the partition types.
   * @param other other SoAView buffer.
   */
  template <size_t... additionalPartitionTypeIndices>
  void appendAdditionalTypesImpl(const SoAViewType &other, std::index_sequence<additionalPartitionTypeIndices...>) {
    // fold expression
    (appendSingleType<additionalPartitionTypeIndices>(other), ...);
  }

  /**
   * Applies append to every partition of a given type with the relevant partition from the other SoA buffer. Increases the maximum
   * depth of the partition type if the other SoA buffer has deeper partitions of this type.
   * @tparam additionalPartitionTypeIndex index corresponding to the desired additional SoA partition's type.
   * @param other other SoA buffer.
   */
  template <size_t additionalPartitionTypeIndex>
  void appendSingleType(SoA<SoAType> &other) {
    // check if resize is needed and if so resize
    const auto otherMaxDepth = other->template getMaxDepthOfGivenType<additionalPartitionTypeIndex>();
    if (otherMaxDepth > this->template getMaxDepthOfGivenType<additionalPartitionTypeIndex>()) {
      this->template resizeMaxDepthOfGivenType<additionalPartitionTypeIndex>(otherMaxDepth);
    }

    // Append SoAs up until otherMaxDepth.
    for (size_t i = 0; i < otherMaxDepth; ++i) {
      std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions)[i].append(std::get<additionalPartitionTypeIndex>(other)[i]);
    }
  }

  /**
   * Applies append to every partition of a given type with the relevant partition from the other SoAView buffer. Increases the maximum
   * depth of the partition type if the other SoAView buffer has deeper partitions of this type.
   * @tparam additionalPartitionTypeIndex index corresponding to the desired additional SoA partition's type.
   * @param other other SoAView buffer.
   */
  template <size_t additionalPartitionTypeIndex>
  void appendSingleType(SoAView<SoAType> &other) {
    using additionalPartitionType = std::tuple_element<additionalPartitionTypeIndex, typename SoAType::AdditionalPartitionTypesReference>;
    // check if resize is needed and if so resize
    const auto otherMaxDepth = other->template getMaxDepthOfGivenType<additionalPartitionTypeIndex>();
    if (otherMaxDepth > this->template getMaxDepthOfGivenType<additionalPartitionTypeIndex>()) {
      this->template resizeMaxDepthOfGivenType<additionalPartitionTypeIndex>(otherMaxDepth);
    }

    // Append SoAs up until otherMaxDepth.
    for (size_t i = 0; i < otherMaxDepth; ++i) {
      // A SoAView is a view on a SoA with normal SoAPartitions. Therefore, the SoAPartitionViews must be constructed from
      // the SoAPartitions.
      const auto partitionView = SoAPartitionView<additionalPartitionType>(std::get<additionalPartitionTypeIndex>(other)[i], other->_s);
      std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions)[i].append(partitionView);
    }
  }

  /**
   * Returns a tuple of pointers to the vectors of the given attributes in the main partition.
   * @tparam attributes parameter pack of the IDs of the desired attributes.
   * @return Tuple of pointers to the beginnings of the desired attribute vectors.
   */
  template <class Functor, size_t... I>
  auto beginMainAttrImpl(std::index_sequence<I...>) {
    return std::make_tuple(_mainSoAPartition.template begin<std::get<I>(Functor::getNeededAttr())>()...);
  }

  /**
   * Implementation of begin for the additional partitions.
   *
   * @tparam additionalAttributesTypes tuple of the types of the returned attribute pointers.
   * @tparam additionalAttr the desired attributes (of type additionalAttributesType)
   * @tparam typeIndices parameter pack of indices representing each SoAPartitionType
   * @return a tuple (an element for each SoAPartitionType) of a vector (of size maxDepth for each SoAPartitionType) of
   * a tuple (an element for each desired attribute of that SoAPartitionType) of pointers to the start of the corresponding
   * arrays.
   */
  template <typename Functor, size_t... additionalPartitionTypeIndices>
  auto beginAdditionalAttrImpl(std::index_sequence<additionalPartitionTypeIndices...>) {
    return std::make_tuple(beginAdditionalAttrOfSingleKind<Functor, additionalPartitionTypeIndices>()...);
  }

  /**
   * Returns a pointer to the start of all arrays of the given desired attributes of one partition type. Pointers are
   * arranged in a vector of tuples of pointers where the size of the vector is the max depth of this type and the
   * size of the tuple is the number of desired attributes of this type.
   *
   * @tparam additionalAttributesType the type of the returned attribute pointers.
   * @tparam additionalAttr the desired attributes (of type additionalAttributesSingleKindType)
   * @tparam additionalPartitionTypeIndex the index of this partition type.
   * @return
   */
  template <typename Functor, size_t additionalPartitionTypeIndex>
  auto beginAdditionalAttrOfSingleKind() {
    using additionalPartitionType = typename std::tuple_element<additionalPartitionTypeIndex, AdditionalPartitionTypesReference>;
    std::vector<additionalPartitionType> &soaPartitionsOfSingleKindPtrs{};

    const auto maxDepthOfGivenType = getMaxDepthOfGivenType<additionalPartitionTypeIndex>();
    soaPartitionsOfSingleKindPtrs.resize(maxDepthOfGivenType);

    const auto noNeededAttributesOfGivenType = std::get<additionalPartitionTypeIndex>(Functor::getNeededAdditionalAttr()).size();

    for (size_t depth = 0; depth < maxDepthOfGivenType; ++depth) {
      auto const pointerPackTmp = beginAdditionalAttrOfSingleKindImpl<Functor, additionalPartitionTypeIndex>(std::make_index_sequence<noNeededAttributesOfGivenType>{});
      for (size_t i = 0; i < std::tuple_size<additionalPartitionType>(); ++i ) {
        soaPartitionsOfSingleKindPtrs[depth][i] = pointerPackTmp[i];
      }
    }
    return soaPartitionsOfSingleKindPtrs;
  }

  template <typename Functor, size_t additionalPartitionTypeIndex, size_t... attr>
  auto beginAdditionalAttrOfSingleKindImpl(size_t depth, std::index_sequence<attr...>) {
    // The partition from which the arrays the returned pointers will point to.
    const auto partitionOfInterest = std::get<additionalPartitionTypeIndex>(_additionalSoAPartitions)[depth];
    // Extract the array of needed attributes that correspond to this type only.
    const auto attributesOfInterest = std::get<additionalPartitionTypeIndex>(Functor::getNeededAdditionalAttr());
    // Return the needed pointers for this partition.
    return std::make_tuple(partitionOfInterest.template begin<std::get<attr>(attributesOfInterest)>()...);
  }

  /**
   * Helper for the actual implementation of getMaxDepths that takes the partition types as indices in a parameter pack.
   * @tparam additionalPartitionTypeIndices parameter pack of indices for each additionalPartitionType
   * @return array of max depths for each additionalPartitionType
   */
  template <size_t... additionalPartitionTypeIndices>
  std::array<size_t, numAdditionalTypes> getMaxDepthsImpl(std::index_sequence<additionalPartitionTypeIndices...>) {
    return std::array<size_t, numAdditionalTypes>{getMaxDepthOfGivenType<additionalPartitionTypeIndices>()...};
  }
};
}  // namespace autopas
