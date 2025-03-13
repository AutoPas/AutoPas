/**
 * @file VCLListIntersectionTraversal3B.h
 * @author mueller
 * @date 07.01.25
 */

#pragma once

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/Cluster.h"
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
      processClusterListIntersection(cluster, false);
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:

  /**
   * Invokes the calculations of all interactions within the cluster and, if they exist, it's neighbors using list
   * intersection
   * @param cluster
   * @param isHaloCluster
   */
  void processClusterListIntersection(internal::Cluster<Particle> &cluster, bool isHaloCluster) {
    if (isHaloCluster) {
      return;
    }
    _clusterFunctor.traverseClusterTriwise(cluster);
    if (*(cluster.getNeighbors()->data())) {
      auto intersectingNeighbors = std::vector<internal::Cluster<Particle>*>();

      auto &neighborClusters1 = *(cluster.getNeighbors());
      std::sort(neighborClusters1.begin(), neighborClusters1.end(), [](const auto cluster1, const auto cluster2)
                {return (*cluster1)[0].getID() < (*cluster2)[0].getID();});
      auto neighborClusters1End = neighborClusters1.end();

      for (auto neighborClusterIter1 = neighborClusters1.begin(); neighborClusterIter1 != neighborClusters1End;) {
        internal::Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
        //take care of all cluster pairs first
        _clusterFunctor.traverseClusterPairTriwise(cluster, neighbor1);

        auto &neighborClusters2 = *(neighbor1.getNeighbors());
        std::sort(neighborClusters2.begin(), neighborClusters2.end(), [](const auto cluster1, const auto cluster2)
                  {return (*cluster1)[0].getID() < (*cluster2)[0].getID();});

        //reserve space in buffer
        std::size_t maxIntersectionSize = std::min(neighborClusters1.size(), neighborClusters2.size());
        intersectingNeighbors.reserve(maxIntersectionSize);

        //preemptively increment neighbor list since current neighbor can not be in intersection
        ++neighborClusterIter1;

        //intersect the neighbors of cluster and neighbor1 to find common neighbors
        std::set_intersection(
            neighborClusterIter1, neighborClusters1End,
            neighborClusters2.begin(), neighborClusters2.end(),
            std::back_inserter(intersectingNeighbors),
            [](const auto cluster1, const auto cluster2){return (*cluster1)[0].getID() < (*cluster2)[0].getID();});

        //iterate over intersected list
        for (auto neighborClusterIter2 = intersectingNeighbors.begin(); neighborClusterIter2 != intersectingNeighbors.end(); ++neighborClusterIter2) {
          internal::Cluster<Particle> &neighbor2 = **(neighborClusterIter2);
          _clusterFunctor.traverseClusterTriplet(cluster, neighbor1, neighbor2);
        }
        //clear buffer
        intersectingNeighbors.clear();
      }
    }
  }

  TriwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, TriwiseFunctor> _clusterFunctor;
};
} // namespace autopas
