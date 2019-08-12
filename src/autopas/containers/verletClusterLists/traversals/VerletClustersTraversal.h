/**
 * @file VerletClustersTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/verletClusterLists/traversals/ClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversalInterface.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.
 * @tparam dataLayout The data layout to use. Currently, only AoS is supported.
 * @tparam useNewton3 If newton 3 should be used. Currently, only false is supported.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class VerletClustersTraversal : public TraversalInterface,
                                public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;
  static constexpr size_t clusterSize = VerletClusterLists<Particle>::clusterSize;

 public:
  /**
   * Constructor of the VerletClustersTraversal.
   * @param pairwiseFunctor The functor to use for the traveral.
   */
  explicit VerletClustersTraversal(PairwiseFunctor *pairwiseFunctor)
      : _functor(pairwiseFunctor), _clusterFunctor(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClusters; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }
  bool getUseNewton3() const override { return useNewton3; }

  bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa) and not useNewton3;
  }

  void initTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    clusterList.loadParticlesIntoSoAs(_functor);
  }

  void endTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    clusterList.extractParticlesFromSoAs(_functor);
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    const auto _clusterTraverseFunctor = [this](internal::Cluster<Particle, clusterSize> &cluster) {
      _clusterFunctor.traverseCluster(cluster);
      for (auto *neighborCluster : cluster.getNeighbors()) {
        _clusterFunctor.traverseClusterPair(cluster, *neighborCluster);
      }
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  PairwiseFunctor *_functor;
  ClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3, clusterSize> _clusterFunctor;
};
}  // namespace autopas
