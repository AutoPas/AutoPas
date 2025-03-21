/**
 * @file OTC01Traversal.h
 *
 * @author Johannes Spies
 * @date 16.06.2021
 */

#pragma once

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/traversals/OTTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/DataLayoutConverter.h"
#include "autopas/utils/logging/OctreeLogger.h"

namespace autopas {

/**
 * This traversal is capable of iterating over particles stored in the Octree data structure.
 *
 * @tparam Particle_T
 * @tparam PairwiseFunctor
 */
template <class Particle_T, class PairwiseFunctor>
class OTC01Traversal : public CellTraversal<OctreeLeafNode<Particle_T>>,
                       public OTTraversalInterface<OctreeNodeWrapper<Particle_T>> {
 public:
  /**
   * A shortcut to specify the type of the actual iterated cell
   */
  using ParticleCell = OctreeLeafNode<Particle_T>;

  /**
   * Constructor for the Octree traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff (this is enough for the octree traversal, please don't use the interaction length here.)
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit OTC01Traversal(PairwiseFunctor *pairwiseFunctor, double cutoff, double interactionLength,
                          DataLayoutOption dataLayout, bool useNewton3)
      : CellTraversal<ParticleCell>({2, 1, 1}),
        OTTraversalInterface<OctreeNodeWrapper<Particle_T>>(interactionLength, dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/, dataLayout,
                     useNewton3),
        _dataLayoutConverter(pairwiseFunctor, dataLayout) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ot_c01; }

  [[nodiscard]] bool isApplicable() const override { return not this->_useNewton3; }

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
   * @copydoc TraversalInterface::traverseParticles()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;

    auto *haloWrapper = this->getHalo();

    OctreeLogger<Particle_T>::octreeToJSON(this->getOwned()->getRaw(), this->getHalo()->getRaw(), this->_ownedLeaves,
                                           this->_haloLeaves);

    // Get neighboring cells for each leaf
    // #pragma omp parallel for
    for (int i = 0; i < this->_ownedLeaves.size(); ++i) {
      OctreeLeafNode<Particle_T> *leaf = this->_ownedLeaves[i];

      // Process cell itself
      _cellFunctor.processCell(*leaf);

      // Process connection to all neighbors in this octree
      auto uniqueNeighboringLeaves = leaf->getNeighborLeaves();
      for (OctreeLeafNode<Particle_T> *neighborLeaf : uniqueNeighboringLeaves) {
        _cellFunctor.processCellPair(*leaf, *neighborLeaf);
      }

      // Process particles in halo cell that are in range
      auto min = leaf->getBoxMin() - this->_interactionLength;
      auto max = leaf->getBoxMax() + this->_interactionLength;
      auto haloNeighbors = haloWrapper->getLeavesInRange(min, max);

      for (OctreeLeafNode<Particle_T> *neighborLeaf : haloNeighbors) {
        _cellFunctor.processCellPair(*leaf, *neighborLeaf);
      }
    }
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellFunctor.setSortingThreshold(sortingThreshold); }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<ParticleCell, PairwiseFunctor, /*bidirectional*/ false> _cellFunctor;

  /**
   * Data Layout Converter to be used with this traversal
   */
  utils::DataLayoutConverter<PairwiseFunctor> _dataLayoutConverter;
};

}  // namespace autopas