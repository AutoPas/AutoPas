/**
 * @file VerletNeighborListInterface.h
 *
 * @date 20.05.2019
 * @author humig
 */

#pragma once

#include <vector>
#include "autopas/options/TraversalOption.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"

namespace autopas {

template<class Particle>
class VerletNeighborListInterface {
  // something with aos to soa

  // TODO: Maybe add some kind of traversal selector
 public:
  virtual ~VerletNeighborListInterface() = default;

  virtual std::vector<TraversalOption> getAllTraversals() = 0;

  virtual void buildNeighborList(LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                                       typename VerletListHelpers<Particle>::SoAArraysType> &linkedCells,
                                 bool useNewton3) = 0;

};

} //namespace autopas