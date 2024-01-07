/**
 * @file VCLSlicedBalancedTraversal.h
 * @author fischerv
 * @date 09 Jun 2020
 */

#pragma once

#include "autopas/containers/cellPairTraversals/SlicedBalancedBasedTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * This traversal splits the domain into slices along the longer dimension among x and y.
 * The slices are processed in parallel by multiple threads. Race conditions are prevented,
 * by placing locks on the starting layers of each slice.
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 */
template <class ParticleCell, class PairwiseFunctor>
class VCLSlicedBalancedTraversal : public SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, false>,
                                   public VCLTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;

  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor> _clusterFunctor;

  void processBaseStep(unsigned long x, unsigned long y) {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
    auto &currentTower = clusterList.getTowerByIndex(x, y);

    for (auto clusterIter =
             (this->_useNewton3 ? currentTower.getClusters().begin() : currentTower.getFirstOwnedCluster());
         clusterIter < (this->_useNewton3 ? currentTower.getClusters().end() : currentTower.getFirstTailHaloCluster());
         ++clusterIter) {
      const auto isHaloCluster =
          clusterIter < currentTower.getFirstOwnedCluster() or clusterIter >= currentTower.getFirstTailHaloCluster();
      _clusterFunctor.processCluster(*clusterIter, isHaloCluster);
    }
  }

 public:
  /**
   * Constructor of the VCLSlicedBalancedTraversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor to use for the traveral.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param clusterSize the number of particles per cluster.
   * @param dataLayout
   * @param useNewton3
   */
  explicit VCLSlicedBalancedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                      const double interactionLength, const std::array<double, 3> &cellLength,
                                      size_t clusterSize, DataLayoutOption::Value dataLayout, bool useNewton3)
      : SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, false>(dims, pairwiseFunctor, interactionLength,
                                                                           cellLength, dataLayout, useNewton3),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize, dataLayout, useNewton3) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_sliced_balanced; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return this->_dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return this->_useNewton3; }

  void loadDataLayout() override {
    if (this->_dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
  }

  void endTraversal() override {
    if (this->_dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
  }

  void traverseParticlePairs() override {
    this->slicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) { processBaseStep(x, y); });
  }

  /**
   * @copydoc autopas::CellPairTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}
};
}  // namespace autopas
