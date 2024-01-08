/**
 * @file VCLClusterIterationTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "autopas/containers/TraversalBase.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.

 */
template <class ParticleCell, class PairwiseFunctor>
class VCLClusterIterationTraversal : public TraversalBase,
                                     public VCLTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor of the VCLClusterIterationTraversal.
   * @param pairwiseFunctor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   * @param dataLayout The data layout to use. Currently, only AoS is supported.
   * @param useNewton3 If newton 3 should be used. Currently, only false is supported.
   */
  explicit VCLClusterIterationTraversal(PairwiseFunctor *pairwiseFunctor, size_t clusterSize,
                                        DataLayoutOption::Value dataLayout, bool useNewton3)
      : TraversalBase(dataLayout, useNewton3),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize, dataLayout, useNewton3) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_cluster_iteration; }

  [[nodiscard]] bool isApplicable() const override {
    return (_dataLayout == DataLayoutOption::aos || _dataLayout == DataLayoutOption::soa) and not _useNewton3;
  }

  void initTraversal() override {
    if (_dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
  }

  void endTraversal() override {
    if (_dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;

    const auto _clusterTraverseFunctor = [this](internal::Cluster<Particle> &cluster) {
      _clusterFunctor.processCluster(cluster, false);
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor> _clusterFunctor;
};
}  // namespace autopas
