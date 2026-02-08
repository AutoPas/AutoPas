/**
 * @file VCLC01BalancedTraversal.h
 * @author humig
 * @date 12.08.19
 */

#pragma once

#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 *
 * It uses a static scheduling that gives each thread about the same amount of cluster pairs to handle.
 * @tparam Particle_T Particle Type.
 * @tparam PairwiseFunctor The type of the functor.
 */
template <class Particle_T, class PairwiseFunctor>
class VCLC01BalancedTraversal : public TraversalInterface, public VCLTraversalInterface<Particle_T> {
 public:
  /**
   * Constructor of the VCLC01BalancedTraversal.
   * @param pairwiseFunctor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   * @param dataLayout The data layout to use.
   * @param useNewton3 If newton 3 should be used. Only false is supported.
   */
  explicit VCLC01BalancedTraversal(PairwiseFunctor *pairwiseFunctor, size_t clusterSize, DataLayoutOption dataLayout,
                                   bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize, dataLayout, useNewton3) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_c01_balanced; }

  [[nodiscard]] bool isApplicable() const override {
    return (_dataLayout == DataLayoutOption::aos or _dataLayout == DataLayoutOption::soa) and not _useNewton3;
  }

  void initTraversal() override {
    if (_dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VCLTraversalInterface<Particle_T>::_verletClusterLists;
    clusterList.loadParticlesIntoSoAs(_functor);
  }

  void endTraversal() override {
    if (_dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VCLTraversalInterface<Particle_T>::_verletClusterLists;
    clusterList.extractParticlesFromSoAs(_functor);
  }

  void traverseParticles() override {
    auto &clusterList = *VCLTraversalInterface<Particle_T>::_verletClusterLists;
    auto &clusterThreadPartition = clusterList.getClusterThreadPartition();

    const auto maxThreads = clusterThreadPartition.size();
    const auto numThreads = std::clamp(maxThreads, 1ul, static_cast<size_t>(autopas_get_preferred_num_threads()));
    AUTOPAS_OPENMP(parallel for num_threads(numThreads))
    for(auto i = 0; i < maxThreads; i++) {
      const auto &clusterRange = clusterThreadPartition[i];
      auto &towers = *VCLTraversalInterface<Particle_T>::_towers;
      size_t clusterCount = 0;
      for (size_t towerIndex = clusterRange.startTowerIndex;
          clusterCount < clusterRange.numClusters and towerIndex < towers.size(); towerIndex++) {
        auto &currentTower = towers[towerIndex];
        auto startIndexInTower =
            clusterCount == 0 ? clusterRange.startIndexInTower : currentTower.getFirstOwnedClusterIndex();
        for (size_t clusterIndex = startIndexInTower;
            clusterCount < clusterRange.numClusters and clusterIndex < currentTower.getFirstTailHaloClusterIndex();
            clusterIndex++, clusterCount++) {
          const auto isHaloCluster = clusterIndex < currentTower.getFirstOwnedClusterIndex() or
                                    clusterIndex >= currentTower.getFirstTailHaloClusterIndex();
          _clusterFunctor.processCluster(currentTower.getCluster(clusterIndex), isHaloCluster);
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
  internal::VCLClusterFunctor<Particle_T, PairwiseFunctor> _clusterFunctor;
};
}  // namespace autopas
