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
#include "autopas/utils/CudaSoA.h"
#include "autopas/utils/SoAView.h"

#if defined(AUTOPAS_CUDA)
#include "autopas/pairwiseFunctors/FunctorCuda.cuh"
#endif

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
   * Get attributes needed for computation.
   * @return Attributes needed for computation.
   * @todo C++20: make this function virtual
   */
  constexpr static std::array<typename Particle::AttributeNames, 0> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 0>{};
  }

  /**
   * Get attributes needed for computation without N3 optimization.
   * @return Attributes needed for computation.
   * @todo C++20: make this function virtual
   */
  constexpr static std::array<typename Particle::AttributeNames, 0> getNeededAttr(std::false_type) {
    return Functor_T::getNeededAttr();
  }

  /**
   * Get attributes computed by this functor.
   * @return Attributes computed by this functor.
   * @todo C++20: make this function virtual
   */
  constexpr static std::array<typename Particle::AttributeNames, 0> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 0>{};
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
                                const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                                bool newton3) {
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
   * Functor using Cuda on SoA in device Memory
   *
   * This Functor calculates the pair-wise interactions between particles in the device_handle on the GPU
   *
   * @param device_handle soa in device memory
   * @param newton3 defines whether or whether not to use newton
   */
  virtual void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle, bool newton3) {
    utils::ExceptionHandler::exception("Functor::CudaFunctorNoN3: not yet implemented");
  }

  /**
   * Functor using Cuda on SoAs in device Memory
   *
   * This Functor calculates the pair-wise interactions between particles in the device_handle1 and device_handle2 on
   * the GPU
   *
   * @param device_handle1 first soa in device memory
   * @param device_handle2 second soa in device memory
   * @param newton3 defines whether or whether not to use newton
   */
  virtual void CudaFunctor(CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle1,
                           CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle2, bool newton3) {
    utils::ExceptionHandler::exception("Functor::CudaFunctorNoN3(two cells): not yet implemented");
  }

  /**
   * Copies the SoA data of the given cell to the Cuda device.
   *
   * @param soa  Structure of arrays where the data is loaded.
   * @param device_handle soa in device memory where the data is copied to
   */
  virtual void deviceSoALoader(SoA<SoAArraysType> &soa,
                               CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) {
    utils::ExceptionHandler::exception("Functor::CudaDeviceSoALoader: not yet implemented");
  }

  /**
   * Copies the data stored on the Cuda device back to the SoA overwrites the data in the soa
   *
   * @param soa  Structure of arrays where the data copied to.
   * @param device_handle soa in device memory where the data is loaded from
   */
  virtual void deviceSoAExtractor(SoA<SoAArraysType> &soa,
                                  CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) {
    utils::ExceptionHandler::exception("Functor::CudaDeviceSoAExtractor: not yet implemented");
  }

  /**
   * Copies the AoS data of the given cell in the given soa.
   *
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added
   * to the SoA with the specified offset.
   * @tparam ParticleCell Type of the cell.
   */
  template <class ParticleCell>
  void SoALoader(ParticleCell &cell, SoA<SoAArraysType> &soa, size_t offset) {
    SoALoaderImpl(cell, soa, offset, std::make_index_sequence<Functor_T::getNeededAttr().size()>{});
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
   * Check whether the given clusterSize is appropriate and can be used by the functor.
   * @param clusterSize The size of the clusters.
   * @param dataLayout The used data layout.
   * @return true, iff the cluster size is appropriate.
   */
  virtual bool isAppropriateClusterSize(unsigned int clusterSize, DataLayoutOption::Value dataLayout) const = 0;

#if defined(AUTOPAS_CUDA)
  /**
   * Provides an interface for traversals to directly access Cuda Functions
   * @return Pointer to CudaWrapper of the Functor
   */
  virtual CudaWrapperInterface<typename Particle::ParticleSoAFloatPrecision> *getCudaWrapper() { return nullptr; }

  /**
   * Creates a Cuda SoA object containing all the relevant pointers from the generic Cuda SoA
   * @param device_handle
   * @return unique pointer to the object
   */
  virtual std::unique_ptr<FunctorCudaSoA<typename Particle::ParticleSoAFloatPrecision>> createFunctorCudaSoA(
      CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) {
    return std::make_unique<FunctorCudaSoA<typename Particle::ParticleSoAFloatPrecision>>();
  }
#endif

  /**
   * Getter for the functor's cutoff
   * @return
   */
  double getCutoff() const { return _cutoff; }

 private:
  /**
   * Implements loading of SoA buffers.
   * @tparam cell_t Cell type.
   * @tparam I Attribute.
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added
   * to the SoA with the specified offset.
   */

  template <typename cell_t, std::size_t... I>
  void SoALoaderImpl(cell_t &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset, std::index_sequence<I...>) {
    soa.resizeArrays(offset + cell.numParticles());

    if (cell.numParticles() == 0) return;

    /**
     * Store the start address of all needed arrays inside the SoA buffer in a tuple. This avoids unnecessary look ups
     * in the following loop.
     */
    // maybe_unused necessary because gcc doesnt understand that pointer is used later
    [[maybe_unused]] auto const pointer = std::make_tuple(soa.template begin<Functor_T::getNeededAttr()[I]>()...);

    auto cellIter = cell.begin();
    // load particles in SoAs
    for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
      /**
       * The following statement writes the values of all attributes defined in neededAttr into the respective position
       * inside the SoA buffer. I represents the index inside neededAttr. The whole expression is folded sizeof...(I)
       * times over the comma operator. E.g. like this (std::index_sequence<I...> = 0, 1):
       * ((std::get<0>(pointer)[i] = cellIter->template get<Functor_T::getNeededAttr()[0]>()),
       * (std::get<1>(pointer)[i] = cellIter->template get<Functor_T::getNeededAttr()[1]>()))
       */
      ((std::get<I>(pointer)[i] = cellIter->template get<Functor_T::getNeededAttr()[I]>()), ...);
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
    if (cell.numParticles() == 0) return;

    /**
     * Store the start address of all needed arrays inside the SoA buffer in a tuple. This avoids unnecessary look ups
     * in the following loop.
     */
    // maybe_unused necessary because gcc doesnt understand that pointer is used later
    [[maybe_unused]] auto const pointer = std::make_tuple(soa.template begin<Functor_T::getComputedAttr()[I]>()...);

    auto cellIter = cell.begin();
    // write values in SoAs back to particles
    for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
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