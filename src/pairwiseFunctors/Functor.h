/*
 * Functor.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_FUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_FUNCTOR_H_

namespace autopas {

template <class Particle>
class Functor {
 public:
  virtual ~Functor() = default;

  /**
   * Functor for arrays of structures (AoS).
   * This functor should calculate the forces or any other pair-wise interaction between two particles.
   * This should include a cutoff check if needed!
   */
  virtual void AoSFunctor(Particle &, Particle &) {}

  /**
   *
   */
  virtual void SoAFunctor() {}
  //	virtual void SoALoader() = 0
  //	virtual void SoAStorer() = 0
  //	virtual void SoAInserter() = 0;
  //	virtual void SoAExtracter() = 0;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_FUNCTOR_H_ */
