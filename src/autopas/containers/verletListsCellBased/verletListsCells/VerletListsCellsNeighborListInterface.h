/**
 * @file VerletListsCellsNeighborListInterface.h
 * @author tirgendetwas
 * @date 27.10.20
 */

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
  /**default destructor*/
    ~VerletListsCellsNeighborListInterface() = default;

  /**
   * Builds AoS neighbor list from underlying linked cells object.
   * @param linkedCells Linked Cells object used to build the neighor list.
   * @param useNewton3 Whether Newton 3 should be used for the neighbor list.
   * @param cutoff
   * @param skin
   * @param interactionLength
   * @param buildTraversalOption Traversal option necessary for generator functor.
   * */
  virtual void buildAoSNeighborList(LinkedCells<typename VerletListsCellsHelpers<Particle>::VLCCellType> &linkedCells, bool useNewton3,
                            double cutoff, double skin, double interactionLength, const TraversalOption buildTraversalOption) = 0;

  /**
   * Get the neighbors list of a particle.
   * @param particle
   * @return the neighbor list of the particle
   */
  virtual const std::vector<Particle *> &getVerletList(const Particle *particle) const = 0;
};

}
