/**
 * @file OTTraversalInterface.h
 *
 * @author Johannes Spies
 * @date 11.05.2021
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/octree/OctreeNodeWrapper.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This interface exists to provide a row interface for octree to add its cells.
 */
template <typename Particle, typename ParticleCell>
class OTTraversalInterface {
 public:
  OTTraversalInterface(double interactionLength) : _interactionLength(interactionLength) {}

  /**
   * Notify the traversal about the cells that it is able to traverse.
   * @param cells A vector of size 2 containing the owned and the halo octrees.
   */
  // virtual void setCells(std::vector<ParticleCell> *cells) = 0;

  /**
   * Set the cells to iterate.
   * @param cells A list of octree roots that should be used during iteration
   */
  void setCells(std::vector<OctreeNodeWrapper<Particle>> *cells) { _cells = cells; }

  long leafGathering, startConversion, endConversion, leafClearing;

 protected:
  template <typename PairwiseFunctor, DataLayoutOption::Value dataLayout>
  void loadBuffers(utils::DataLayoutConverter<PairwiseFunctor, dataLayout> &dataLayoutConverter,
                   OctreeNodeWrapper<Particle> *wrapper, std::vector<OctreeLeafNode<Particle> *> &leaves) {
    TIME_IT(leafGathering, wrapper->appendAllLeaves(leaves));

    TIME_IT(startConversion, for (OctreeLeafNode<Particle> *leaf : leaves) {
      dataLayoutConverter.loadDataLayout(*leaf);
    });
  }

  template <typename PairwiseFunctor, DataLayoutOption::Value dataLayout>
  void unloadBuffers(utils::DataLayoutConverter<PairwiseFunctor, dataLayout> &dataLayoutConverter,
                     std::vector<OctreeLeafNode<Particle> *> &leaves) {
    TIME_IT(endConversion, for (OctreeLeafNode<Particle> *leaf : leaves) {
      dataLayoutConverter.storeDataLayout(*leaf);
    });

    // Remove the cached leaves
    TIME_IT(leafClearing, leaves.clear());
  }

  OctreeNodeWrapper<Particle> *getOwned() { return dynamic_cast<OctreeNodeWrapper<Particle> *>(&(*_cells)[0]); }

  OctreeNodeWrapper<Particle> *getHalo() { return dynamic_cast<OctreeNodeWrapper<Particle> *>(&(*_cells)[1]); }

  std::vector<OctreeNodeWrapper<Particle>> *_cells;

  /**
   * A list of all leaves in the owned octree
   */
  std::vector<OctreeLeafNode<Particle> *> _ownedLeaves;

  /**
   * A list of all leaves in the halo octree
   */
  std::vector<OctreeLeafNode<Particle> *> _haloLeaves;

  /**
   * The interaction length is used for finding neighbors
   */
  double _interactionLength;
};

}  // namespace autopas
