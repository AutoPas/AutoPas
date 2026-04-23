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
template <class Particle_T, class CRTP_T>
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

  /**
   * SoAFunctorPair with interaction sorting along a projection axis.
   * Implementations project all particles onto sortingDirection, sort both views by projection (ascending),
   * then iterate with early break when |proj_j - proj_i| > sortingCutoff.
   *
   * @param soa1 First SoA view.
   * @param soa2 Second SoA view.
   * @param sortingDirection Normalized vector along the cell-pair axis.
   * @param sortingCutoff 1D cutoff for early break (normally the interaction cutoff).
   * @param newton3 Whether to apply Newton's third law.
   */
  virtual void SoAFunctorPairSorted(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2,
                                    const std::array<double, 3> & /*sortingDirection*/, double /*sortingCutoff*/,
                                    bool newton3) {
    SoAFunctorPair(soa1, soa2, newton3);
  }

  /**
   * Specifies whether the functor is capable of using the specified Vectorization Pattern in the SoA functor.
   *
   * Note: All pairwise functors should support p1xVec by default.
   * @param vecPattern
   * @return whether the functor is capable of using the specified Vectorization Pattern
   */
  bool isVecPatternAllowed(const VectorizationPatternOption::Value vecPattern) override {
    return vecPattern == VectorizationPatternOption::p1xVec;
  }
};

}  // namespace autopas