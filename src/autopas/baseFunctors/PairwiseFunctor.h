/**
 * @file PairwiseFunctor.h
 *
 * @date 12.08.2023
 * @author muehlhaeusser
 */

#pragma once

#include <type_traits>

#include "Functor.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/KokkosStorage.h"
#include "autopas/utils/SoAView.h"

namespace autopas {

template <class Particle>
class VerletListHelpers;

/**
 * PairwiseFunctor class. This class describes the pairwise interactions between
 * particles.
 * @copydoc autopas::Functor
 *
 * @tparam Particle_T the type of Particle
 * @tparam CRTP_T the actual type of the functor
 */
template <class Particle_T, class CRTP_T, class MemSpace = Kokkos::HostSpace>
class PairwiseFunctor : public Functor<Particle_T, CRTP_T> {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Constructor
   * @param cutoff
   */
  explicit PairwiseFunctor(double cutoff) : Functor<Particle_T, CRTP_T>(cutoff){};

  virtual ~PairwiseFunctor() = default;

  /**
   * PairwiseFunctor for arrays of structures (AoS).
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between two particles.
   * This should include a cutoff check if needed!
   * @param i Particle i
   * @param j Particle j
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) {
    utils::ExceptionHandler::exception("{}::AoSFunctor: not implemented", this->getName());
  }

  virtual void AoSFunctorKokkos(Particle_T &i, Particle_T &j, bool newton3) {
    // TODO: implement
  }

  /**
   * PairwiseFunctor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles in an soa.
   * This should include a cutoff check if needed!
   *
   * @param soa Structure of arrays
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) {
    utils::ExceptionHandler::exception("{}::SoAFunctorSingle: not implemented", this->getName());
  }

  KOKKOS_INLINE_FUNCTION
  virtual void SoAFunctorSingleKokkos(int i, const Particle_T::KokkosSoAArraysType& soa, size_t N, bool newton3) {
    // No Op unless overridden
  }

  /**
   * PairwiseFunctor for structure of arrays (SoA) for neighbor lists
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
    utils::ExceptionHandler::exception("{}::SoAFunctorVerlet: not implemented", this->getName());
  }

  /**
   * PairwiseFunctor for structure of arrays (SoA)
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
    utils::ExceptionHandler::exception("{}::SoAFunctorPair: not implemented", this->getName());
  }

  KOKKOS_INLINE_FUNCTION
  virtual void SoAFunctorPairKokkos(Particle_T::KokkosSoAArraysType& soa1, Particle_T::KokkosSoAArraysType& soa2, bool newton3) {
    // No Op unless overridden
  }
};

}  // namespace autopas