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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VCLSlicedBalancedTraversal
    : public SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
      public VCLTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;

  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3> _clusterFunctor;

  void processBaseStep(unsigned long x, unsigned long y) {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
    auto &currentTower = clusterList.getTowerAtCoordinates(x, y);
    for (auto &cluster : currentTower.getClusters()) {
      _clusterFunctor.traverseCluster(cluster);
      for (auto *neighborCluster : cluster.getNeighbors()) {
        _clusterFunctor.traverseClusterPair(cluster, *neighborCluster);
      }
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
   */
  explicit VCLSlicedBalancedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                      const double interactionLength, const std::array<double, 3> &cellLength,
                                      size_t clusterSize)
      : SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(
            dims, pairwiseFunctor, interactionLength, cellLength),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_sliced_balanced; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  void initTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
    SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::initTraversal();
  }

  void endTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
    SlicedBalancedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::endTraversal();
  }

  void traverseParticlePairs() override {
    this->template slicedTraversal</*allCells*/ true>(
        [&](unsigned long x, unsigned long y, unsigned long z) { processBaseStep(x, y); });
  }
};
}  // namespace autopas
