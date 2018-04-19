/*
 * VerletLists.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_VERLETLISTS_H_
#define SRC_CONTAINERS_VERLETLISTS_H_

#include "LinkedCells.h"

namespace autopas {

template <class Particle, class ParticleCell>
class VerletLists : public LinkedCells<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the VerletLists class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   */
  VerletLists(const std::array<double, 3> boxMin,
              const std::array<double, 3> boxMax, double cutoff)
      : LinkedCells<Particle, ParticleCell>(boxMin, boxMax, cutoff) {}

 private:
  // ThreeDimensionalCellHandler
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_VERLETLISTS_H_ */
