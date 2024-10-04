/**
 * @file Functor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <type_traits>

#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/SoAView.h"
#include "autopas/utils/logging/FLOPLogger.h"

namespace autopas {

/**
 * Newton 3 modes for the LJFunctor.
 */
enum class FunctorN3Modes {
  Newton3Only,
  Newton3Off,
  Both,
};

template <class Particle>
class VerletListHelpers;

//template <class Particle, class CRTP_T>
//class Functor;

/**
 * Functor class. This class describes the pairwise interactions between
 * particles.
 * Both an array of structure (AoS) and a structure of array (SoA) are supported
 * to be used with functors.
 * Newton3: A functor does not have to implement both a newton3 and a
 * non-newton3 version. Instead you can specify, which version you use by
 * overriding allowsNonNewton3 resp. allowsNewton3
 *
 * @tparam Particle the type of Particle
 * @tparam ParticleCell_t the type of ParticleCell
 */
template <class Particle, class CRTP_T>
class Functor {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Make the Implementation type template publicly available.
   */
  using Functor_T = CRTP_T;

  /**
   * Constructor
   * @param cutoff
   */
  explicit Functor(double cutoff) : _cutoff(cutoff){};

  virtual ~Functor() = default;

  /**
   * This function is called at the start of each traversal.
   * Use it for resetting global values or initializing them.
   */
  virtual void initTraversal(){};

  /**
   * This function is called at the end of each traversal.
   * You may accumulate values in this step.
   * @param newton3
   */
  virtual void endTraversal(bool newton3){};

  /**
   * Functor for arrays of structures (AoS).
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between two particles.
   * This should include a cutoff check if needed!
   * @param i Particle i
   * @param j Particle j
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void AoSFunctor(Particle &i, Particle &j, bool newton3) {
    utils::ExceptionHandler::exception("Functor::AoSFunctor: not yet implemented");
  }


  /**
   * Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles in an soa.
   * This should include a cutoff check if needed!
   *
   * @param soa Structure of arrays
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) {
    utils::ExceptionHandler::exception("Functor::SoAFunctorSingle: not yet implemented");
  }

  /**
   * Functor for structure of arrays (SoA) for neighbor lists
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between the particle in the SoA with index indexFirst and all particles with indices in the neighborList.
   * This should include a cutoff check if needed!
   *
   * @param soa Structure of arrays
   * @param indexFirst The index of the first particle for each interaction
   * @param neighborList The list of neighbors
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                                const std::vector<size_t, AlignedAllocator<size_t>> &neighborList, bool newton3) {
    utils::ExceptionHandler::exception("Functor::SoAFunctorVerlet: not yet implemented");
  }

  /**
   * Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles of soa1 and soa2.
   * This should include a cutoff check if needed!
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3) {
    utils::ExceptionHandler::exception("Functor::SoAFunctorPair: not yet implemented");
  }

  /**
   * Copies the AoS data of the given cell in the given soa.
   *
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added to the SoA with the specified offset.
   * @param skipSoAResize If resizing of the SoA buffers should be skipped or not. If this is called with true, it must
   * be ensured before the call that there is sufficient capacity in the SoA.
   * @note The parameter skipSoAResize is usually set to false, only for VerletListsCellBased Containers it is set to
   * true, since they resize the SoA before the call to SoALoader.
   * @tparam ParticleCell Type of the cell.
   */
  template <class ParticleCell>
  void SoALoader(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset, bool skipSoAResize) {
    using getNeededAttrN3ReturnType = typename decltype(std::function{Functor_T::getNeededAttr()})::result_type;

    using getComputedAttrReturnType = typename decltype(std::function{Functor_T::getComputedAttr()})::result_type;

    using getNeededAdditionalAttrReturnType = typename decltype(std::function{Functor_T::getNeededAdditionalAttr()})::result_type;

    static constexpr auto numAdditionalPartitions = Functor_T::getNeededAdditionalAttr().size();

    if (not skipSoAResize) {
      soa.resizeArrays(offset + cell.size());
    }

    if (cell.isEmpty()) return;

    /**
     * Store the start address of all needed arrays inside the SoA buffer in a tuple. This avoids unnecessary look ups
     * in the following loop.
     */
    // maybe_unused necessary because gcc doesn't understand that pointer is used later
    [[maybe_unused]] const auto [mainPtr, additionalPtr] = soa.template begin<getNeededAttrN3ReturnType , Functor_T::getNeededAttr(),
                                                    getNeededAdditionalAttrReturnType, Functor_T::getNeededAdditionalAttr()>();

    // get number of partitions of each type
    const auto maxDepthsOfAdditionalPartitionTypes = soa.getMaxDepths();

    auto cellIter = cell.begin();

    // some aliases that are required as template parameters. These are simply the types of mainPtr, additionalPtr, and
    // cellIter i.e. the result_types of the functions that are used to obtain these pointers.
    using mainPtrType = typename decltype(std::function{
        std::get<0>(soa.template begin<getNeededAttrN3ReturnType , Functor_T::getNeededAttr(),
                                       getNeededAdditionalAttrReturnType, Functor_T::getNeededAdditionalAttr()>())})::result_type;
    using additionalPtrType = typename decltype(std::function{
        std::get<1>(soa.template begin<getNeededAttrN3ReturnType , Functor_T::getNeededAttr(),
                                       getNeededAdditionalAttrReturnType, Functor_T::getNeededAdditionalAttr()>())})::result_type;
    using cellIterType = typename decltype(std::function{cell.begin()})::result_type;

    // load particles one-by-one into the SoAPartitions
    for (size_t i = offset; cellIter != cell.end(); ++cellIter, ++i) {
      SoALoaderMainPartitionImpl<cellIterType, mainPtrType>(cellIter, mainPtr, i, std::make_index_sequence<Functor_T::getNeededAttr().size()>{});
      SoALoaderAdditionalPartitionsImpl<cellIterType, additionalPtrType>(cellIter, additionalPtr, i, maxDepthsOfAdditionalPartitionTypes, std::make_index_sequence<numAdditionalPartitions>{});
    }
  }

  /**
   * Copies the data stored in the soa back into the cell.
   *
   * @param cell Cell where the data should be stored.
   * @param soa  Structure of arrays from where the data is loaded.
   * @param offset Offset within the SoA. The data of the soa should be
   * extracted starting at offset.
   * @tparam ParticleCell Type of the cell.
   */
  template <typename ParticleCell>
  void SoAExtractor(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset) {
    SoAExtractorImpl(cell, soa, offset, std::make_index_sequence<Functor_T::getComputedAttr().size()>{});
  }

  /**
   * Specifies whether the functor is capable of Newton3-like functors.
   * If the functor provides an interface to soa or aos functions that utilize
   * Newton's third law of motion (actio = reactio) to reduce the computational
   * complexity this function should return true. If this is not the case this
   * function should return false.
   * @return true if and only if this functor provides an interface to
   * Newton3-like functions.
   */
  virtual bool allowsNewton3() = 0;

  /**
   * Specifies whether the functor is capable of non-Newton3-like functors.
   * If the functor provides an interface to soa or aos functions that do not
   * utilize Newton's third law of motion (actio = reactio) this function should
   * return true. If this is not the case this function should return false.
   * @return true if and only if this functor provides an interface to functions
   * that do not utilize Newton3.
   */
  virtual bool allowsNonNewton3() = 0;

  /**
   * Specifies whether the functor should be considered for the auto-tuning process.
   * @return true if and only if this functor is relevant for auto-tuning.
   */
  virtual bool isRelevantForTuning() = 0;

  /**
   * Getter for the functor's cutoff
   * @return
   */
  double getCutoff() const { return _cutoff; }

  /**
   * Get the number of FLOPs. Implementation required if FLOPLogger used.
   *
   * If derived class provides no implementation, the FLOPLogger interprets the default numeric_limits<size_t>::max()
   * output as invalid and leaves "Not Implemented" the log.
   * @return number of FLOPs
   */
  [[nodiscard]] virtual size_t getNumFLOPs() const { return std::numeric_limits<size_t>::max(); }

  /**
   * Get the hit rate. Implementation required if FLOPLogger used.
   *
   * If derived class provides no implementation, the FLOPLogger interprets the default NaN output as invalid and
   * leaves "Not Implemented" in the log.
   *
   * @return (number of kernel calls) / (number of distance calculations)
   */
  [[nodiscard]] virtual double getHitRate() const { return std::numeric_limits<double>::quiet_NaN(); }

 private:
  /**
   * Implementation of SoALoader for the main partition. Uses a fold expression to "loop" over each attribute. Loads
   * only a single particle.
   * @tparam cellIterType type of cellIter
   * @tparam mainPtrType type of mainPtr. Should be a tuple of pointers to arrays of the relevant types.
   * @tparam attr parameter pack of  0,1,2,... for each desired attributes.
   * @param cellIter iterator pointing to particle in the cell.
   * @param mainPtr pointer structure of type mainPtrType.
   * @param arrayIndex index of particle (that is being loaded) in each array of the SoA.
   */
  template <typename cellIterType, typename mainPtrType, size_t... attr>
  void SoALoaderMainPartitionImpl(cellIterType &cellIter, mainPtrType mainPtr, size_t arrayIndex, std::index_sequence<attr...>) {
    /**
     * The following statement writes the values of all attributes defined in neededAttr into the respective position
     * inside the SoA buffer. attr represents the index inside neededAttr. The whole expression is folded sizeof...(attr)
     * times over the comma operator. E.g. like this (std::index_sequence<attr...> = 0, 1):
     * ((std::get<0>(mainPtr)[arrayIndex] = cellIter->template get<Functor_T::getNeededAttr()[0]>()),
     * (std::get<1>(mainPtr)[arrayIndex] = cellIter->template get<Functor_T::getNeededAttr()[1]>()))
     */
    ((std::get<attr>(mainPtr)[arrayIndex] = cellIter->template get<Functor_T::getNeededAttr()[attr]>()), ...);
  }

  /**
   * Implementation of SoALoader for the additional partitions. Uses a fold expression to "loop" over each additional
   * partition type. Loads only a single particle.
   * @tparam cellIterType type of cellIter
   * @tparam additionalPtrType type of additionalPtr. Should be a tuple of a vector of a tuple of pointers to arrays of
   * the relevant types.
   * @tparam additionalPartitionTypeIndices parameter pack of indices of additional partition types.
   * @param cellIter iterator pointing to particle in the cell.
   * @param additionalPtr pointer structure of type additionalPtrType.
   * @param arrayIndex index of particle (that is being loaded) in each array of the SoA.
   * @param maxDepthsOfAdditionalPartitionTypes array of maximum depths of each additional partition type (<=> the size
   * of the vector of SoAPartitions of that type)
   */
  template <typename cellIterType, typename additionalPtrType, size_t... additionalPartitionTypeIndices>
  void SoALoaderAdditionalPartitionsImpl(cellIterType &cellIter, additionalPtrType additionalPtr, size_t arrayIndex,
                                         std::array<size_t, sizeof...(additionalPartitionTypeIndices)> maxDepthsOfAdditionalPartitionTypes,
                                         std::index_sequence<additionalPartitionTypeIndices...>) {
    (SoALoaderAdditionalPartitionsSingleType<cellIterType, additionalPtrType, additionalPartitionTypeIndices>(
         cellIter, additionalPtr, arrayIndex, maxDepthsOfAdditionalPartitionTypes[additionalPartitionTypeIndices]),...);
  }

  /**
   * Loads all SoAPartitions of a single additional partition type. Loads only a single particle.
   * @tparam cellIterType type of cellIter
   * @tparam additionalPtrType type of additionalPtr. Should be a tuple of a vector of a tuple of pointers to arrays of
   * the relevant types.
   * @tparam additionalPartitionTypeIndex index of the additional partition type.
   * @param cellIter iterator pointing to particle in the cell.
   * @param additionalPtr pointer structure of type additionalPtrType.
   * @param arrayIndex index of particle (that is being loaded) in each array of the SoA.
   * @param maxDepth maximum depth of this additional partition type.
   */
  template <typename cellIterType, typename additionalPtrType, size_t additionalPartitionTypeIndex>
  void SoALoaderAdditionalPartitionsSingleType(cellIterType &cellIter, additionalPtrType additionalPtr, size_t arrayIndex,
                                               size_t maxDepth) {
    constexpr size_t numNeededAdditionalAttr = std::get<additionalPartitionTypeIndex>(Functor_T::getNeededAdditionalAttr()).size();

    SoALoaderAdditionalPartitionsSingleTypeImpl<cellIterType, additionalPtrType, additionalPartitionTypeIndex>(
        cellIter, additionalPtr, arrayIndex, maxDepth, std::make_index_sequence<numNeededAdditionalAttr>{});
  }

  /**
   * Implementation of SoALoaderAdditionalPartitionsSingleType. Uses a fold expression to "loop" over each attribute of
   * the additional partition type. Uses a conventional for-loop to loop over each partition of that additional partition
   * type. Loads only a single particle.
   * @tparam cellIterType type of cellIter
   * @tparam additionalPtrType type of additionalPtr. Should be a tuple of a vector of a tuple of pointers to arrays of
   * the relevant types.
   * @tparam additionalPartitionTypeIndex index of the additional partition type.
   * @tparam attr parameter pack of indices 0,1,2,... for each attribute being loaded.
   * @param cellIter iterator pointing to particle in the cell.
   * @param additionalPtr pointer structure of type additionalPtrType.
   * @param arrayIndex index of particle (that is being loaded) in each array of the SoA.
   * @param maxDepth maximum depth of this additional partition type.
   */
  template <typename cellIterType, typename additionalPtrType, size_t additionalPartitionTypeIndex, size_t... attr>
  void SoALoaderAdditionalPartitionsSingleTypeImpl(cellIterType &cellIter, additionalPtrType additionalPtr,
                                                   size_t arrayIndex, size_t maxDepth, std::index_sequence<attr...>) {
    // Loop over each partition of this type
    for (size_t depth = 0; depth < maxDepth; ++depth) {
      /**
       * This is a fold expression that "loop" over each attribute (see SoALoaderMainPartitionImpl). For each attribute,
       * it copies the RHS to the LHS:
       * * RHS: Calls the particle's getter with desired additionalPartitionTypeIndex and desired additional attribute.
       * * LHS: This is the "arrayIndex"'th element of the "attr"'th array in the relevant SoAPartition. This SoAPartition
       *   has type index "additionalPartitionTypeIndex" and is the partition of depth "depth" of this type.
       */
      ((std::get<attr>(std::get<additionalPartitionTypeIndex>(additionalPtr)[depth])[arrayIndex] =
            cellIter->template get<std::get<additionalPartitionTypeIndex, (Functor_T::getNeededAdditionalAttr())[attr]>(depth)>), ...);
    }

  }

  /**
   * Implements extraction of SoA buffers.
   * @tparam cell_t Cell type.
   * @tparam I Attribute.
   * @param cell Cell.
   * @param soa SoA buffer.
   * @param offset Offset
   */
  template <typename cell_t, std::size_t... I>
  void SoAExtractorImpl(cell_t &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset, std::index_sequence<I...>) {
    if (cell.isEmpty()) return;

    /**
     * Store the start address of all needed arrays inside the SoA buffer in a tuple. This avoids unnecessary look ups
     * in the following loop.
     */
    // maybe_unused necessary because gcc doesn't understand that pointer is used later
    [[maybe_unused]] auto const pointer = std::make_tuple(soa.template begin<Functor_T::getComputedAttr()[I]>()...);

    auto cellIter = cell.begin();
    // write values in SoAs back to particles
    for (size_t i = offset; cellIter != cell.end(); ++cellIter, ++i) {
      /**
       * The following statement writes the value of all attributes defined in computedAttr back into the particle.
       * I represents the index inside computedAttr.
       * The whole expression is folded sizeof...(I) times over the comma operator. E.g. like this
       * (std::index_sequence<I...> = 0, 1):
       * (cellIter->template set<Functor_T::getComputedAttr()[0]>(std::get<0>(pointer)[i]),
       * cellIter->template set<Functor_T::getComputedAttr()[1]>(std::get<1>(pointer)[i]))
       */
      (cellIter->template set<Functor_T::getComputedAttr()[I]>(std::get<I>(pointer)[i]), ...);
    }
  }

  double _cutoff;
};

}  // namespace autopas