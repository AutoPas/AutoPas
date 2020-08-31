/**
 * @file VCLC01BalancedTraversal.h
 * @author humig
 * @date 12.08.19
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 *
 * It uses a static scheduling that gives each thread about the same amount of cluster pairs to handle.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.
 * @tparam dataLayout The data layout to use.
 * @tparam useNewton3 If newton 3 should be used. Only false is supported.
 */
template <class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VCLC01BalancedTraversal : public TraversalInterface, public VCLTraversalInterface<Particle> {
 public:
  /**
   * Constructor of the VCLC01BalancedTraversal.
   * @param pairwiseFunctor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   */
  explicit VCLC01BalancedTraversal(PairwiseFunctor *pairwiseFunctor, size_t clusterSize)
      : _functor(pairwiseFunctor), _clusterFunctor(pairwiseFunctor, clusterSize) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_c01_balanced; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos or dataLayout == DataLayoutOption::soa) and not useNewton3;
  }

  void initTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
    clusterList.loadParticlesIntoSoAs(_functor);
  }

  void endTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
    clusterList.extractParticlesFromSoAs(_functor);
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
    auto &clusterThreadPartition = clusterList.getClusterThreadPartition();

    auto numThreads = clusterThreadPartition.size();
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel num_threads(numThreads)
#endif
    {
      auto threadNum = autopas_get_thread_num();
      const auto &clusterRange = clusterThreadPartition[threadNum];
      auto &towers = *VCLTraversalInterface<Particle>::_towers;
      size_t clusterCount = 0;
      for (size_t towerIndex = clusterRange.startTowerIndex;
           clusterCount < clusterRange.numClusters and towerIndex < towers.size(); towerIndex++) {
        auto &currentTower = towers[towerIndex];
        auto startIndexInTower = clusterCount == 0 ? clusterRange.startIndexInTower : 0;
        for (size_t clusterIndex = startIndexInTower;
             clusterIndex < currentTower.getNumClusters() && clusterCount < clusterRange.numClusters;
             clusterIndex++, clusterCount++) {
          auto &currentCluster = currentTower.getCluster(clusterIndex);
          _clusterFunctor.traverseCluster(currentCluster);
          for (auto *neighborCluster : currentCluster.getNeighbors()) {
            _clusterFunctor.traverseClusterPair(currentCluster, *neighborCluster);
          }
        }
      }
      if (clusterCount != clusterRange.numClusters) {
        autopas::utils::ExceptionHandler::exception(
            "VCLC01BalancedTraversal::traverseParticlePairs(): Not all or too many clusters traversed, probably "
            "the clusterThreadPartitions are wrong! TraversedClusters={}, ClustersInRange={}",
            clusterCount, clusterRange.numClusters);
      }
    }
  }

  bool needsStaticClusterThreadPartition() override { return true; }

 private:
  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3> _clusterFunctor;
};
}  // namespace autopas
