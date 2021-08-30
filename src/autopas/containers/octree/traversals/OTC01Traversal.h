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
#include "autopas/utils/logging/OctreeLogger.h"

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
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   */
  explicit OTC01Traversal(PairwiseFunctor *pairwiseFunctor, double cutoff, double interactionLength)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        OTTraversalInterface<Particle, OctreeNodeWrapper<Particle>>(interactionLength),
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
        _dataLayoutConverter(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ot_c01; }

  [[nodiscard]] bool isApplicable() const override { return not useNewton3; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    // Preprocess all leaves
    this->loadBuffers(_dataLayoutConverter, this->getOwned(), this->_ownedLeaves);
    this->loadBuffers(_dataLayoutConverter, this->getHalo(), this->_haloLeaves);
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
    //OctreeLogger<Particle>::octreeToJSON(this->getOwned()->getRaw(), this->getHalo()->getRaw(), this->_ownedLeaves,
    //                                     this->_haloLeaves);

    // Get neighboring cells for each leaf
    //#pragma omp parallel for default(none) shared(haloWrapper)
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < this->_ownedLeaves.size(); ++i) {
      OctreeLeafNode<Particle> *leaf = this->_ownedLeaves[i];

      // Process cell itself
      _cellFunctor.processCell(*leaf);

      // Process connection to all neighbors in this octree
      auto uniqueNeighboringLeaves = leaf->getNeighborLeaves();
      for (OctreeLeafNode<Particle> *neighborLeaf : uniqueNeighboringLeaves) {
        _cellFunctor.processCellPair(*leaf, *neighborLeaf);
      }

      // Process particles in halo cell that are in range
      auto min = utils::ArrayMath::subScalar(leaf->getBoxMin(), this->_interactionLength);
      auto max = utils::ArrayMath::addScalar(leaf->getBoxMax(), this->_interactionLength);
      auto haloNeighbors = this->getHalo()->getLeavesInRange(min, max);

      for (OctreeLeafNode<Particle> *neighborLeaf : haloNeighbors) {
        _cellFunctor.processCellPair(*leaf, *neighborLeaf);
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