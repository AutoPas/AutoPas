/**
 * @file VCLTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.
 * @tparam dataLayout The data layout to use. Currently, only AoS is supported.
 * @tparam useNewton3 If newton 3 should be used. Currently, only false is supported.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VCLTraversal : public TraversalInterface, public VCLTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor of the VCLTraversal.
   * @param pairwiseFunctor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   */
  explicit VCLTraversal(PairwiseFunctor *pairwiseFunctor, size_t clusterSize)
      : _functor(pairwiseFunctor), _clusterFunctor(pairwiseFunctor, clusterSize) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_cluster_iteration; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa) and not useNewton3;
  }

  void initTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
  }

  void endTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;

    const auto _clusterTraverseFunctor = [this](internal::Cluster<Particle> &cluster) {
      _clusterFunctor.traverseCluster(cluster);
      for (auto *neighborCluster : cluster.getNeighbors()) {
        _clusterFunctor.traverseClusterPair(cluster, *neighborCluster);
      }
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3> _clusterFunctor;
};
}  // namespace autopas
