/**
 * @file VCLListIntersectionTraversal3B.h
 * @author mueller
 * @date 07.01.25
 */

#pragma once

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 * @tparam ParticleCell
 * @tparam TriwiseFunctor The type of the functor.
*/

template <class ParticleCell, class TriwiseFunctor>
class VCLListIntersectionTraversal3B : public TraversalInterface,
                                     public VCLTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor of the VCLListIntersectionTraversal3B.
   * @param functor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   * @param dataLayout The data layout to use. Currently, only AoS is supported.
   * @param useNewton3 If newton 3 should be used. Currently, only false is supported.
   */
  explicit VCLListIntersectionTraversal3B(TriwiseFunctor *functor, size_t clusterSize,
                                        DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        _functor(functor),
        _clusterFunctor(functor, clusterSize, dataLayout, useNewton3) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_list_intersection_3b; }

  [[nodiscard]] bool isApplicable() const override {
    return (_dataLayout == DataLayoutOption::aos or _dataLayout == DataLayoutOption::soa) and not _useNewton3;
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

  void traverseParticles() override {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;

    const auto _clusterTraverseFunctor = [this](internal::Cluster<Particle> &cluster) {
      _clusterFunctor.processClusterListIntersection(cluster, false);
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  TriwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, TriwiseFunctor> _clusterFunctor;
};
} // namespace autopas
