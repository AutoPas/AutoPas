/*
 * Functor.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_FUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_FUNCTOR_H_

#include <utils/SoA.h>

namespace autopas {

/**
 * Functor class. This class describes the pairwise interactions between
 * particles.
 * Both an array of structure (AoS) and a structure of array (SoA) are supported
 * to be used with functors.
 * @tparam Particle the type of Particle
 * @tparam ParticleCell the type of ParticleCell
 */
template <class Particle, class ParticleCell>
class Functor {
 public:
  virtual ~Functor() = default;

  /**
   * @brief Functor for arrays of structures (AoS).
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between two particles.
   * This should include a cutoff check if needed!
   */
  virtual void AoSFunctor(Particle &, Particle &) {}

  /**
   * @brief Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles in an soa.
   * This should include a cutoff check if needed!
   *
   * @param soa Structure of arrays
   */
  virtual void SoAFunctor(SoA &soa) {}

  /**
   * @brief Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction
   * between all particles of soa1 and soa2.
   * This should include a cutoff check if needed!
   *
   * @param soa1 First structure of arrays.
   * @param soa2 Second structure of arrays.
   */
  virtual void SoAFunctor(SoA &soa1, SoA &soa2) {}

  /**
   * @brief Copies the AoS data of the given cell in the given soa.
   *
   * @param cell Cell from where the data is loaded.
   * @param soa  Structure of arrays where the data is copied to.
   */
  virtual void SoALoader(ParticleCell &cell, SoA *soa) {}

  /**
   * @brief Copies the data stored in the soa in the cell.
   *
   * @param cell Cell where the data should be stored.
   * @param soa  Structure of arrays from where the data is loaded.
   */
  virtual void SoAExtractor(ParticleCell *cell, SoA *soa) {}
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_FUNCTOR_H_ */
