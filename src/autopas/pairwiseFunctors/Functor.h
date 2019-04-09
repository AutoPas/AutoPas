/**
 * @file Functor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/CudaSoA.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#endif

namespace autopas {

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
 * @tparam ParticleCell the type of ParticleCell
 */
template <class Particle, class ParticleCell, class SoAArraysType = typename Particle::SoAArraysType>
class Functor {
 public:
  virtual ~Functor() = default;

  /**
   * @brief Functor for arrays of structures (AoS).
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
   * @brief Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles in an soa.
   * This should include a cutoff check if needed!
   *
   * @param soa Structure of arrays
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) {
    utils::ExceptionHandler::exception("Functor::SoAFunctor(one soa): not yet implemented");
  }

  /**
   * @brief Functor for structure of arrays (SoA) for neighbor lists
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between the particles in the SoA that are marked by the verlet list
   * This should include a cutoff check if needed!
   *
   * iFrom an iTo define the range inside of the neighborList that should be
   * iterated over. The starting index is i = iFrom. The iteration will continue
   * while i < iTo.
   *
   * @param soa Structure of arrays
   * @param neighborList The list of neighbors
   * @param iFrom the starting index of the vector neighborList that should be
   * iterated over
   * @param iTo the first index that should not be iterated over. (Should be at
   * least iFrom and less than soa.size())
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctor(SoA<SoAArraysType> &soa,
                          const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList,
                          size_t iFrom, size_t iTo, bool newton3) {
    utils::ExceptionHandler::exception("Functor::SoAFunctor(verlet): not yet implemented");
  }

  /**
   * @brief Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles of soa1 and soa2.
   * This should include a cutoff check if needed!
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3) {
    utils::ExceptionHandler::exception("Functor::SoAFunctor(two soa): not yet implemented");
  }

  /**
   * @brief Functor using Cuda on SoA in device Memory
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
   * @brief Functor using Cuda on SoAs in device Memory
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
   * @brief Copies the SoA data of the given cell to the Cuda device.
   *
   * @param soa  Structure of arrays where the data is loaded.
   * @param device_handle soa in device memory where the data is copied to
   */
  virtual void deviceSoALoader(::autopas::SoA<SoAArraysType> &soa,
                               CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) {
    utils::ExceptionHandler::exception("Functor::CudaDeviceSoALoader: not yet implemented");
  }

  /**
   * @brief Copies the data stored on the Cuda device back to the SoA overwrites the data in the soa
   *
   * @param soa  Structure of arrays where the data copied to.
   * @param device_handle soa in device memory where the data is loaded from
   */
  virtual void deviceSoAExtractor(::autopas::SoA<SoAArraysType> &soa,
                                  CudaSoA<typename Particle::CudaDeviceArraysType> &device_handle) {
    utils::ExceptionHandler::exception("Functor::CudaDeviceSoAExtractor: not yet implemented");
  }

  /**
   * @brief Copies the AoS data of the given cell in the given soa.
   *
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   * @param offset Offset within the SoA. The data of the cell should be added
   * to the SoA with the specified offset.
   */
  virtual void SoALoader(ParticleCell &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset = 0) {
    utils::ExceptionHandler::exception("Functor::SoALoader: not yet implemented");
  }

  /**
   * @brief Copies the data stored in the soa back into the cell.
   *
   * @param cell Cell where the data should be stored.
   * @param soa  Structure of arrays from where the data is loaded.
   * @param offset Offset within the SoA. The data of the soa should be
   * extracted starting at offset.
   */
  virtual void SoAExtractor(ParticleCell &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset = 0) {
    utils::ExceptionHandler::exception("Functor::SoAExtractor: not yet implemented");
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
  virtual bool allowsNewton3() { return true; }

  /**
   * Specifies whether the functor is capable of non-Newton3-like functors.
   * If the functor provides an interface to soa or aos functions that do not
   * utilize Newton's third law of motion (actio = reactio) this function should
   * return true. If this is not the case this function should return false.
   * @return true if and only if this functor provides an interface to functions
   * that do not utilize Newton3.
   */
  virtual bool allowsNonNewton3() { return false; }

  /**
   * Specifies whether the functor should be considered for the auto-tuning process.
   * @return true if and only if this functor is relevant for auto-tuning.
   */
  virtual bool isRelevantForTuning() = 0;

#if defined(AUTOPAS_CUDA)
  /**
   * Provides an interface for traversals to directly access Cuda Functions
   * @return Pointer to CudaWrapper of the Functor
   */
  virtual CudaWrapper *getCudaWrapper() { return NULL; }
#endif
};

/**
 * Macro to define the SoALoaders. The body to be executed is the variadic argument.
 * @param cell name for cell like parameter
 * @param soa name for soa
 * @param offset name for offset
 *
 * @note The last Argument is variadic such that commas pose no problem.
 * @note generates two loaders, one for verlet lists, one for the normal case.
 * @note the need for this could be removed if the soa's are removed from the particlecells (highly unlikely)
 */
#define AUTOPAS_FUNCTOR_SOALOADER(cell, soa, offset, ...)                                                            \
  void SoALoader(ParticleCell &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset = 0) override { __VA_ARGS__ } \
  /** @copydoc SoALoader(ParticleCell &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset) */                   \
  template <typename = std::enable_if_t<not std::is_same<                                                            \
                typename ::autopas::VerletListHelpers<Particle>::VerletListParticleCellType, ParticleCell>::value>>  \
  void SoALoader(typename ::autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,                  \
                 ::autopas::SoA<SoAArraysType> &soa, size_t offset = 0) {                                            \
    __VA_ARGS__                                                                                                      \
  }

/**
 * Macro to define the SoAExtractors. The body to be executed is the variadic argument.
 * @param cell name for cell like parameter
 * @param soa name for soa
 * @param offset name for offset
 *
 * @note ParticleCell and SoAArraysType need defined via typedefs in the functor.
 * @note The last Argument is variadic such that commas pose no problem.
 * @note generates two extractors, one for verlet lists, one for the normal case.
 * @note the need for this could be removed if the soa's are removed from the particlecells (highly unlikely)
 */
#define AUTOPAS_FUNCTOR_SOAEXTRACTOR(cell, soa, offset, ...)                                                        \
  void SoAExtractor(ParticleCell &cell, ::autopas::SoA<SoAArraysType> &soa, size_t offset = 0) override {           \
    __VA_ARGS__                                                                                                     \
  }                                                                                                                 \
  /** @copydoc SoAExtractor(ParticleCell &, ::autopas::SoA<SoAArraysType> &, size_t) */                             \
  template <typename = std::enable_if_t<not std::is_same<                                                           \
                typename ::autopas::VerletListHelpers<Particle>::VerletListParticleCellType, ParticleCell>::value>> \
  void SoAExtractor(typename ::autopas::VerletListHelpers<Particle>::VerletListParticleCellType &cell,              \
                    ::autopas::SoA<SoAArraysType> &soa, size_t offset = 0) {                                        \
    __VA_ARGS__                                                                                                     \
  }

}  // namespace autopas

#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
