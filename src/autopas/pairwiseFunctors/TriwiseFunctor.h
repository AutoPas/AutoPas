/**
 * @file TriwiseFunctor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <type_traits>

#include "Functor.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/SoAView.h"

namespace autopas {

/**
 * TriwiseFunctor class. This class describes the pairwise interactions between
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
class TriwiseFunctor : public Functor<Particle, CRTP_T> {
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
  explicit TriwiseFunctor(double cutoff) : Functor<Particle, CRTP_T>(cutoff){};

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
  virtual void AoSFunctor(Particle &i, Particle &j, Particle &k, bool newton3) {
    utils::ExceptionHandler::exception("TriwiseFunctor::AoSFunctor: not yet implemented");
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
    utils::ExceptionHandler::exception("TriwiseFunctor::SoAFunctorSingle: not yet implemented");
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
    utils::ExceptionHandler::exception("TriwiseFunctor::SoAFunctorVerlet: not yet implemented");
  }

  /**
   * TriwiseFunctor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other triwise interaction
   * between all particles of soa1 and soa2.
   * This should include a cutoff check if needed!
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   * @param newton3 defines whether or whether not to use newton 3
   */
  virtual void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3) {
    utils::ExceptionHandler::exception("TriwiseFunctor::SoAFunctorPair: not yet implemented");
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
  virtual void SoAFunctorTriple(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, SoAView<SoAArraysType> soa3, bool newton3) {
    utils::ExceptionHandler::exception("TriwiseFunctor::SoAFunctorTriple: not yet implemented");
  }

  /**
   * Returns name of functor. Intended for use with the iteration logger, to differentiate between calls to computeInteractions
   * using different functors in the logs.
   * @return name of functor.
   */
  virtual std::string getName() { return "TriwiseFunctor"; }

  /**
   * Return number of interacting bodies. Required to determine the relevant traversals.
   * @return number of interacting bodies of the functor
   */
  static constexpr unsigned int getNBody() { return 3; }
};

}  // namespace autopas