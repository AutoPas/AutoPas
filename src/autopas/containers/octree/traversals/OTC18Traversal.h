/**
 * @file OTNaiveTraversal.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/traversals/OTTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This traversal is capable of iterating over particles stored in the Octree data structure. This traversal does not
 * use any parallelization or speed-increasing strategies and is therefore called naive.
 *
 * @tparam Particle
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class OTC18Traversal : public CellPairTraversal<OctreeLeafNode<Particle>>,
                       public OTTraversalInterface<Particle, OctreeNodeWrapper<Particle>> {
 public:
  /**
   * A shortcut to specify the type of the actual iterated cell
   */
  using ParticleCell = OctreeLeafNode<Particle>;

  /**
   * Constructor for the Octree traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the octree traversal, please don't use the interaction length here.)
   */
  explicit OTC18Traversal(PairwiseFunctor *pairwiseFunctor, double cutoff, double interactionLength)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        OTTraversalInterface<Particle, OctreeNodeWrapper<Particle>>(interactionLength),
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
        _dataLayoutConverter(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ot_naive; }

  [[nodiscard]] bool isApplicable() const override {
    return useNewton3;
    // @todo Re-enable this traversal, when fixing https://github.com/AutoPas/AutoPas/issues/621
    return false;
  }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    // Preprocess all leaves
    this->loadBuffers(_dataLayoutConverter, this->getOwned(), this->_ownedLeaves);
    this->loadBuffers(_dataLayoutConverter, this->getHalo(), this->_haloLeaves);

    // Assign IDs to the owned leaves
    for (int i = 0; i < this->_ownedLeaves.size(); ++i) {
      this->_ownedLeaves[i]->setID(i);
    }

    // Assign IDs to the halo leaves
    for (int i = 0; i < this->_haloLeaves.size(); ++i) {
      this->_haloLeaves[i]->setID(this->_ownedLeaves.size() + i);
    }
  }

  void endTraversal() override {
    // Postprocess all leaves
    this->unloadBuffers(_dataLayoutConverter, this->_ownedLeaves);
    this->unloadBuffers(_dataLayoutConverter, this->_haloLeaves);
  }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override {
    auto *haloWrapper = this->getHalo();

    // Get neighboring cells for each leaf
    for (OctreeLeafNode<Particle> *leaf : this->_ownedLeaves) {
      // Process cell itself
      _cellFunctor.processCell(*leaf);

      // Process connection to all neighbors
      auto uniqueNeighboringLeaves = leaf->getNeighborLeaves();
      for (OctreeLeafNode<Particle> *neighborLeaf : uniqueNeighboringLeaves) {
        if (leaf->getID() < neighborLeaf->getID()) {
          // Execute the cell functor
          _cellFunctor.processCellPair(*leaf, *neighborLeaf);
        }
      }

      // Process particles in halo cell that are in range
      auto min = utils::ArrayMath::subScalar(leaf->getBoxMin(), this->_interactionLength);
      auto max = utils::ArrayMath::addScalar(leaf->getBoxMax(), this->_interactionLength);
      auto haloNeighbors = haloWrapper->getLeavesInRange(min, max);

      for (OctreeLeafNode<Particle> *neighborLeaf : haloNeighbors) {
        if (leaf->getID() < neighborLeaf->getID()) {
          _cellFunctor.processCellPair(*leaf, *neighborLeaf);
        }
      }
    }
  }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<Particle, ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false> _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;
};

}  // namespace autopas