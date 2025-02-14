/**
 * @file VCLTripletListIterationTraversal3B.h
 * @author mueller
 * @date 07.01.25
 */

#pragma once

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas{

/**
 * Traversal for VerletClusterLists. Does not support newton 3.
 * @tparam ParticleCell
 * @tparam TriwiseFunctor The type of the functor.

*/

template <class ParticleCell, class TriwiseFunctor>
class VCLPairListIterationTraversal3B : public TraversalInterface, public VCLTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor of the VCLPairListIterationTraversal3B.
   * @param functor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   * @param dataLayout The data layout to use. Currently, only AoS is supported.
   * @param useNewton3 If newton 3 should be used. Currently, only false is supported.
   */
  explicit VCLPairListIterationTraversal3B(TriwiseFunctor *functor, size_t clusterSize,
                                          DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        _functor(functor),
        _clusterFunctor(functor, clusterSize, dataLayout, useNewton3) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_pair_list_iteration_3b; }

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
      processClusterPairListInteration(cluster, false);
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  /**
   * Invokes the calculations of all interactions within the cluster and, if they exist, it's neighbors using pair lists
   * @param cluster
   * @param isHaloCluster
   */
  void processClusterPairListInteration(internal::Cluster<Particle> &cluster, bool isHaloCluster) {
    //TODO: switch depending on data layout and N3
    if (not isHaloCluster) {
      _clusterFunctor.traverseClusterTriwise(cluster);
    }
    if (*(cluster.getNeighbors()->data())) {
      auto &neighborClusters = *(cluster.getNeighbors());
      auto neighborClustersEnd = neighborClusters.end();
      for (auto neighborClusterIter1 = neighborClusters.begin(); neighborClusterIter1 != neighborClustersEnd; ++neighborClusterIter1) {
        internal::Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
        // 2 options involving 2 clusters:
        // - 1 particle in cluster + 2 particles in neighbor cluster
        // - 2 particles in cluster + 1 particle in neighbor cluster
        _clusterFunctor.traverseClusterPairTriwise(cluster, neighbor1);
      }
    }
    if (cluster.getNeighborPairs()->data()) {
      //traverse pair neighbor list for interactions involving 3 clusters
      auto &neighborClusterPairs = *(cluster.getNeighborPairs());
      auto neighborClusterPairsEnd = neighborClusterPairs.end();
      for (auto neighborClusterPairIter = neighborClusterPairs.begin();
           neighborClusterPairIter != neighborClusterPairsEnd; ++neighborClusterPairIter) {
        internal::Cluster<Particle> &neighbor1 = *((*neighborClusterPairIter).first);
        internal::Cluster<Particle> &neighbor2 = *((*neighborClusterPairIter).second);
        // 1 option involving all 3 clusters:
        // - one particle in cluster, neighbor cluster 1 and neighbor cluster 2 each
        _clusterFunctor.traverseClusterTriplet(cluster, neighbor1, neighbor2);
      }
    }

  }

  TriwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, TriwiseFunctor> _clusterFunctor;
};

} // namespace autopas