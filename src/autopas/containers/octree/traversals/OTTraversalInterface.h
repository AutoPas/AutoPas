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
template <typename ParticleCell>
class OTTraversalInterface {
  using ParticleType = typename ParticleCell::ParticleType;

 public:
  OTTraversalInterface(double interactionLength) : _interactionLength(interactionLength) {}

  /**
   * Notify the traversal about the cells that it is able to traverse.
   * @param cells A vector of size 2 containing the owned and the halo octrees.
   */
  void setCells(std::vector<OctreeNodeWrapper<ParticleType>> *cells) { _cells = cells; }

 protected:
  /**
   * Gather all leaves and load the SoA/AoS buffers.
   *
   * @param dataLayoutConverter The converter to convert the buffers
   * @param wrapper The octree to load the leaves from
   * @param leaves The list to store the leaves in
   */
  template <typename PairwiseFunctor, DataLayoutOption::Value dataLayout>
  void loadBuffers(utils::DataLayoutConverter<PairwiseFunctor, dataLayout> &dataLayoutConverter,
                   OctreeNodeWrapper<ParticleType> *wrapper, std::vector<OctreeLeafNode<ParticleType> *> &leaves) {
    wrapper->appendAllLeaves(leaves);

    for (OctreeLeafNode<ParticleType> *leaf : leaves) {
      dataLayoutConverter.loadDataLayout(*leaf);
    }
  }

  /**
   * Unload the SoA/AoS buffers and clear the gathered leaves list.
   *
   * @param dataLayoutConverter The converter to convert the buffers
   * @param leaves The list to unload the leaves from
   */
  template <typename PairwiseFunctor, DataLayoutOption::Value dataLayout>
  void unloadBuffers(utils::DataLayoutConverter<PairwiseFunctor, dataLayout> &dataLayoutConverter,
                     std::vector<OctreeLeafNode<ParticleType> *> &leaves) {
    for (OctreeLeafNode<ParticleType> *leaf : leaves) {
      dataLayoutConverter.storeDataLayout(*leaf);
    }

    // Remove the cached leaves
    leaves.clear();
  }

  OctreeNodeWrapper<ParticleType> *getOwned() { return dynamic_cast<OctreeNodeWrapper<ParticleType> *>(&(*_cells)[0]); }

  OctreeNodeWrapper<ParticleType> *getHalo() { return dynamic_cast<OctreeNodeWrapper<ParticleType> *>(&(*_cells)[1]); }

  std::vector<OctreeNodeWrapper<ParticleType>> *_cells;

  /**
   * A list of all leaves in the owned octree
   */
  std::vector<OctreeLeafNode<ParticleType> *> _ownedLeaves;

  /**
   * A list of all leaves in the halo octree
   */
  std::vector<OctreeLeafNode<ParticleType> *> _haloLeaves;

  /**
   * The interaction length is used for finding neighbors
   */
  double _interactionLength;
};

}  // namespace autopas
