/**
 * @file VerletNeighborListInterface.h
 *
 * @date 20.05.2019
 * @author humig
 */

#pragma once

#include <vector>
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/TraversalOption.h"

namespace autopas {

template <class Particle>
class VerletNeighborListInterface {
 public:
  virtual ~VerletNeighborListInterface() = default;

  virtual ContainerOption getContainerType() const = 0;

  virtual void buildNeighborList(LinkedCells<Particle, typename VerletListHelpers<Particle>::VerletListParticleCellType,
                                             typename VerletListHelpers<Particle>::SoAArraysType> &linkedCells,
                                 bool useNewton3) = 0;

  virtual void generateSoAFromAoS() = 0;

  virtual bool isSoAListValid() const = 0;

  virtual long getNumberOfNeighborPairs() const = 0;
};

}  // namespace autopas