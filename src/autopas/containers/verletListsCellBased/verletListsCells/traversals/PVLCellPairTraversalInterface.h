/**
 * @file VLCCellPairTraversalInterface.h
 * @author tirgendetwas
 * @date 11.01.21
 */

#pragma once
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/PVLCellPairNeighborList.h"
namespace autopas {
template <class Particle>
/**
 * Interface for traversals used with VLCCellPairNeighborList. Allows the distinction of traversals that are only usable
 * for VLCCellPairNeighborList and not compatible with the VLCAllCellsNeighborList within the VerletListsCells
 * container.
 * @tparam Particle type of particle
 */
class PVLCellPairTraversalInterface {
 public:
  /**
   * Sets the verlet list for the traversal to iterate over.
   * @param verlet The verlet list to iterate over.
   */
  void setPseudoVerletList(PVLCellPairNeighborList<Particle> &pseudoVerlet) { _cellPairPseudoVerletList = &pseudoVerlet; }

 protected:
  /**
   * The verlet list to iterate over.
   */
  PVLCellPairNeighborList<Particle> *_cellPairPseudoVerletList;
};
}  // namespace autopas
