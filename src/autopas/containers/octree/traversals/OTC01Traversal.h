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
                       public OTTraversalInterface<OctreeNodeWrapper<Particle>> {
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
        _cellFunctor(pairwiseFunctor, cutoff /*should use cutoff here, if not used to build verlet-lists*/),
        _dataLayoutConverter(pairwiseFunctor),
        _interactionLength(interactionLength) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::ot_c01; }

  [[nodiscard]] bool isApplicable() const override { return not useNewton3; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; };

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; };

  void initTraversal() override {
    // Preprocess all leaves
    loadBuffers(getOwned(), _ownedLeaves);
    loadBuffers(getHalo(), _haloLeaves);
  }

  void endTraversal() override {
    // Postprocess all leaves
    unloadBuffers(_ownedLeaves);
    unloadBuffers(_haloLeaves);
  }

  /**
   * @copydoc TraversalInterface::traverseParticlePairs()
   * @note This function expects a vector of exactly two cells. First cell is the main region, second is halo.
   */
  void traverseParticlePairs() override {
    auto *haloWrapper = getHalo();

#if 0
    // FOR DEBUGGING ONLY
    // Log all owned leaves for this octree
    fclose(OctreeLogger::leavesToJSON(fopen("owned.json", "w"), _ownedLeaves));
    // Log all halo leaves for this octree
    fclose(OctreeLogger::leavesToJSON(fopen("halo.json", "w"), _haloLeaves));
    FILE *particles = fopen("particles.json", "w");
    fprintf(particles, "{");
    OctreeLogger::particlesToJSON(particles, "owned", getOwned()->getRaw());
    fprintf(particles, ",\n");
    OctreeLogger::particlesToJSON(particles, "halo", getHalo()->getRaw());
    fprintf(particles, "}");
    fclose(particles);
#endif

    // Get neighboring cells for each leaf
    //#pragma omp parallel for
    for (int i = 0; i < _ownedLeaves.size(); ++i) {
      OctreeLeafNode<Particle> *leaf = _ownedLeaves[i];

      // Process cell itself
      _cellFunctor.processCell(*leaf);

      // Process connection to all neighbors in this octree
      auto uniqueNeighboringLeaves = leaf->getNeighborLeaves();
      for (OctreeLeafNode<Particle> *neighborLeaf : uniqueNeighboringLeaves) {
        _cellFunctor.processCellPair(*leaf, *neighborLeaf);
      }

      // Process particles in halo cell that are in range
      auto min = utils::ArrayMath::subScalar(leaf->getBoxMin(), _interactionLength);
      auto max = utils::ArrayMath::addScalar(leaf->getBoxMax(), _interactionLength);
      auto haloNeighbors = haloWrapper->getLeavesInRange(min, max);

      for (OctreeLeafNode<Particle> *neighborLeaf : haloNeighbors) {
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
  void loadBuffers(OctreeNodeWrapper<Particle> *wrapper, std::vector<OctreeLeafNode<Particle> *> &leaves) {
    wrapper->appendAllLeaves(leaves);

    for (OctreeLeafNode<Particle> *leaf : leaves) {
      leaf->clearAlreadyProcessedList();
      _dataLayoutConverter.loadDataLayout(*leaf);
    }
  }

  void unloadBuffers(std::vector<OctreeLeafNode<Particle> *> &leaves) {
    for (OctreeLeafNode<Particle> *leaf : leaves) {
      _dataLayoutConverter.storeDataLayout(*leaf);
    }

    // Remove the cached leaves
    leaves.clear();
  }

  OctreeNodeWrapper<Particle> *getOwned() { return dynamic_cast<OctreeNodeWrapper<Particle> *>(&(*_cells)[0]); }

  OctreeNodeWrapper<Particle> *getHalo() { return dynamic_cast<OctreeNodeWrapper<Particle> *>(&(*_cells)[1]); }

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