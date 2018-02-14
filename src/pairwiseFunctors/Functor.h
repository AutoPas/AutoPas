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

template <class Particle>
class Functor {
 public:
  virtual ~Functor() = default;

  /**
   * @brief Functor for arrays of structures (AoS).
   *
   * This functor should calculate the forces or any other pair-wise interaction between two particles.
   * This should include a cutoff check if needed!
   */
  virtual void AoSFunctor(Particle &, Particle &) {}

  /**
   * @brief Functor for structure of arrays (SoA)
   *
   * This functor should calculate the forces or any other pair-wise interaction between all particles of soa1 and soa2.
   * This should include a cutoff check if needed!
   *
   * @param soa1 first structure of arrays
   * @param soa2 second structure of arrays
   */
  virtual void SoAFunctor(SoA &soa1, SoA &soa2) {}
  // TODO: in what form should particles be passed to this?
  virtual void SoALoader(std::vector<Particle> &particles, SoA *soa) {};
  //	virtual void SoAStorer() = 0
  //	virtual void SoAInserter() = 0;
  virtual void SoAExtractor(std::vector<Particle> *particles, SoA *soa) {};
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_FUNCTOR_H_ */
