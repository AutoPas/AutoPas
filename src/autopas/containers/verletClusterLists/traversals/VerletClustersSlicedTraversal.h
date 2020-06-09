/**
 * @file VerletClustersSlicedTraversal.h
 * @author fischerv
 * @date 09 Jun 2020
 */

#pragma once

#include "autopas/containers/cellPairTraversals/SlicedBasedTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/ClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversalInterface.h"

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
class VerletClustersSlicedTraversal
    : public SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
      public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;

  PairwiseFunctor *_functor;
  internal::ClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3, VerletClusterLists<Particle>::clusterSize>
      _clusterFunctor;

  void processBaseStep(unsigned long x, unsigned long y) {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
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
   * Constructor of the VerletClustersSlicedTraversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor to use for the traveral.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit VerletClustersSlicedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                         const double interactionLength, const std::array<double, 3> &cellLength)
      : SlicedBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                    interactionLength, cellLength),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::verletClustersSliced; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  void initTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VerletClustersTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
  }

  void endTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VerletClustersTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
  }

  void traverseParticlePairs() override {
    this->template slicedTraversal</*allCells*/ true>(
        [&](unsigned long x, unsigned long y, unsigned long z) { processBaseStep(x, y); });
  }
};
}  // namespace autopas
