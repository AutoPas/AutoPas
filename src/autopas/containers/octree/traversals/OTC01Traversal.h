/**
 * @file OTC01Traversal.h
 *
 * @author Johannes Spies
 * @date 16.06.2021
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
 * This traversal is capable of iterating over particles stored in the Octree data structure.
 *
 * @tparam Particle
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class OTC01Traversal : public CellPairTraversal<OctreeLeafNode<Particle>>,
                       public OTTraversalInterface<OctreeNodeWrapper<Particle>> {
 public:
  /**
   * A shortcut to specify the type of the actual iterated cell
   */
  using ParticleCell = OctreeLeafNode<Particle>;

  // TODO(johannes): The TraversalSelector passes the interactionLength as the cutoff value: Keep in mind when
  // implementing...
  /**
   * Constructor for the Octree traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the octree traversal, please don't use the interaction length here.)
   */
  explicit OTC01Traversal(PairwiseFunctor *pairwiseFunctor, double cutoff)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
        _dataLayoutConverter(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ot_c01; }

  [[nodiscard]] bool isApplicable() const override { return not useNewton3; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    // Gather all leaves
    auto *wrapper = dynamic_cast<OctreeNodeWrapper<Particle> *>(&(*_cells)[0]);
    wrapper->appendAllLeaves(_leaves);

    // Preprocess all leaves
    for (OctreeLeafNode<Particle> *leaf : _leaves) {
      leaf->clearAlreadyProcessedList();
      _dataLayoutConverter.loadDataLayout(*leaf);
    }
  }

  void endTraversal() override {
    // Postprocess all leaves
    for (OctreeLeafNode<Particle> *leaf : _leaves) {
      _dataLayoutConverter.storeDataLayout(*leaf);
    }

    // Remove the cached leaves
    _leaves.clear();
  }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override {
    // Get neighboring cells for each leaf
    for (OctreeLeafNode<Particle> *leaf : _leaves) {
      // Process cell itself
      _cellFunctor.processCell(*leaf);

      // Process connection to all neighbors
      auto uniqueNeighboringLeaves = leaf->getNeighborLeaves();
      for (OctreeLeafNode<Particle> *neighborLeaf : uniqueNeighboringLeaves) {
        _cellFunctor.processCellPair(*leaf, *neighborLeaf);
      }
    }
  }

  /**
   * Set the cells to iterate.
   * @param cells A list of octree roots that should be used during iteration
   */
  void setCells(std::vector<OctreeNodeWrapper<Particle>> *cells) override { _cells = cells; }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<Particle, ParticleCell, PairwiseFunctor, dataLayout, useNewton3, false> _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor, dataLayout> _dataLayoutConverter;

  std::vector<OctreeNodeWrapper<Particle>> *_cells;

  /**
   * A list of all leaves in the octree
   * TODO: Include the halo octree
   */
  std::vector<OctreeLeafNode<Particle> *> _leaves;
};

}  // namespace autopas