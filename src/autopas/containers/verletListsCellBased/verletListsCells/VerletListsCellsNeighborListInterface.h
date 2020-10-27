//
// Created by TinaVl on 10/27/2020.
//

#pragma once

#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/containers/linkedCells/LinkedCells.h"

namespace autopas
{

template<class Particle>
class VerletListsCellsNeighborListInterface
{
 public:
  /**TODO*/
    ~VerletListsCellsNeighborListInterface() = default;

  /**TODO*/
  virtual void buildAoSNeighborList(LinkedCells<typename VerletListsCellsHelpers<Particle>::VLCCellType> &linkedCells, bool useNewton3,
                            double cutoff, double skin, double interactionLength, const TraversalOption buildTraversalOption) = 0;

  /**TODO*/
  virtual const std::vector<Particle *> &getVerletList(const Particle *particle) const = 0;
};

}
