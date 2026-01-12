/**
 * @file TriwiseFunctor.h
 *
 * @date 12.08.2023
 * @author muehlhaeusser
 */

#pragma once

#include <type_traits>

#include "Functor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/SoAView.h"

namespace autopas {

/**
 * TriwiseFunctor class. This class describes the triwise interactions between
 * particles.
 * @copydoc autopas::Functor
 *
 * @tparam Particle_T the type of Particle
 * @tparam CRTP_T the actual type of the functor
 */
template <class Particle_T, class CRTP_T>
class TriwiseFunctor : public Functor<Particle_T, CRTP_T> {
 public:
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Constructor
   * @param cutoff
   */
  explicit TriwiseFunctor(double cutoff) : Functor<Particle_T, CRTP_T>(cutoff){};

  virtual ~TriwiseFunctor() = default;

  /**
   * TriwiseFunctor for arrays of structures (AoS).
   *
   * This functor should calculate the forces or any triwise interaction
   * between three particles.
   * This should include a cutoff check if needed!
   * @param i Particle i
   * @param j Particle j
   * @param k Particle k
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void AoSFunctor(Particle_T &i, Particle_T &j, Particle_T &k, bool newton3) {
    utils::ExceptionHandler::exception("{}::AoSFunctor: not implemented", this->getName());
  }

  /**
   * TriwiseFunctor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other triwise interaction
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
   * TriwiseFunctor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other triwise interaction
   * between all particles of soa1 and soa2. It should always calculate forces for all particles in soa1, even when
   * newton3 == false.
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
   * TriwiseFunctor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other triwise interaction
   * between all particles of soa1 and soa2 and soa3.
   * This should include a cutoff check if needed!
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param soa3 Third structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctorTriple(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, SoAView<SoAArraysType> soa3,
                                bool newton3) {
    utils::ExceptionHandler::exception("{}::SoAFunctorTriple: not implemented", this->getName());
  }

  /**
   * TriwiseFunctor for structure of arrays (SoA) for neighbor lists
   *
   * This functor should calculate the forces or any other triwise interaction
   * between the particles in the SoA with index indexFirst and all particles with indices in the neighborList.
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
   * Specifies whether the functor is capable of using the specified Vectorization Pattern in the SoA functor.
   *
   * Note: Currently Vectorization Patterns are not implemented for threebody interactions. p1xVec is used as default.
   * @param vecPattern
   * @return whether the functor is capable of using the specified Vectorization Pattern
   */
  bool isVecPatternAllowed(const VectorizationPatternOption::Value vecPattern) override {
    return vecPattern == VectorizationPatternOption::p1xVec;
  }
};

}  // namespace autopas