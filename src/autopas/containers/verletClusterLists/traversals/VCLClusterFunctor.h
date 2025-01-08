/**
 * @file VCLClusterFunctor.h
 * @author humig, mueller
 * @date 12.08.19
 */

#pragma once

#include "autopas/containers/verletClusterLists/Cluster.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/checkFunctorType.h"

namespace autopas::internal {

/**
 * Provides methods to traverse a single cluster and a pair of clusters.
 *
 * @tparam Particle The type of particle the clusters contain.
 * @tparam Functor The type of the functor the VCLClusterFunctor should use.
 */
template <class Particle, class Functor>
class VCLClusterFunctor {
 public:
  /**
   * Constructs a VCLClusterFunctor that uses the given functor internally.
   * @param functor The functor to use internally.
   * @param clusterSize Number of particles per cluster.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit VCLClusterFunctor(Functor *functor, size_t clusterSize, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(functor), _clusterSize(clusterSize), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Invokes the calculations of all interactions within the cluster and, if they exist, it's neighbors.
   * @param cluster
   * @param isHaloCluster
   */
  void processCluster(internal::Cluster<Particle> &cluster, bool isHaloCluster) {
    if constexpr (utils::isPairwiseFunctor<Functor>()) {
      if (not isHaloCluster) {
        traverseClusterPairwise(cluster);
      }
      // only iterate neighbors if the neighbor list contains more than just nullptr
      if (*(cluster.getNeighbors()->data())) {
        for (auto *neighborClusterPtr : *(cluster.getNeighbors())) {
          traverseClusterPair(cluster, *neighborClusterPtr);
        }
      }
    } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
      if (not isHaloCluster) {
        traverseClusterTriwise(cluster);
      }
      if (*(cluster.getNeighbors()->data())) {
        auto &neighborClusters = *(cluster.getNeighbors());
        auto neighborClustersEnd = neighborClusters.end();
        for (auto neighborClusterIter1 = neighborClusters.begin(); neighborClusterIter1 != neighborClustersEnd; ++neighborClusterIter1) {
          for (auto neighborClusterIter2 = neighborClusterIter1 + 1; neighborClusterIter2 != neighborClustersEnd; ++neighborClusterIter2) {
            Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
            Cluster<Particle> &neighbor2 = **(neighborClusterIter2);
            traverseClusterTriplet(cluster, neighbor1, neighbor2);
          }
        }
      }
    } else {
      utils::ExceptionHandler::exception(
          "VCLClusterFunctor::processCluster(): Functor {} is not of type PairwiseFunctor or TriwiseFunctor.",
          _functor->getName());
    }
  }

  /**
   * Invokes the calculations of all interactions within the cluster and, if they exist, it's neighbors using list
   * intersection
   * @param cluster
   * @param isHaloCluster
   */
  void processClusterListIntersection(internal::Cluster<Particle> &cluster, bool isHaloCluster) {
    if (not isHaloCluster) {
      traverseClusterTriwise(cluster);
    }
    if (*(cluster.getNeighbors()->data())) {
      auto &neighborClusters = *(cluster.getNeighbors());
      //sorting first neighbor list necessary?
      auto neighborClustersEnd = neighborClusters.end();
      for (auto neighborClusterIter1 = neighborClusters.begin(); neighborClusterIter1 != neighborClustersEnd; ++neighborClusterIter1) {
        Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
        auto &neighborClusters1 = *(neighbor1.getNeighbors());
        //sorting second neighbor list necessary?
        //intersect them
        auto intersectedNeighborsIter = neighborClusterIter1;
        auto intersectedNeighborsEnd = std::set_intersection(++intersectedNeighborsIter, neighborClustersEnd, neighborClusters1.begin(), neighborClusters1.end());
        //iterate over intersected list
        for (auto neighborClusterIter2 = intersectedNeighborsIter; neighborClusterIter2 != intersectedNeighborsEnd; ++neighborClusterIter2) {
          Cluster<Particle> &neighbor2 = **(neighborClusterIter1);
          traverseClusterTriplet(cluster, neighbor1, neighbor2);
        }
      }
    }
  }

 private:
  /**
   * Traverses pairs of all particles in the given cluster. Always uses newton 3 in the AoS data layout.
   * @param cluster The cluster to traverse.
   */
  void traverseClusterPairwise(internal::Cluster<Particle> &cluster) {
    if (_dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = i + 1; j < _clusterSize; j++) {
          // this if else branch is needed because of https://github.com/AutoPas/AutoPas/issues/426
          if (_useNewton3) {
            _functor->AoSFunctor(cluster[i], cluster[j], true);
          } else {
            _functor->AoSFunctor(cluster[i], cluster[j], false);
            _functor->AoSFunctor(cluster[j], cluster[i], false);
          }
        }
      }
    } else {
      _functor->SoAFunctorSingle(cluster.getSoAView(), _useNewton3);
    }
  }

  /**
   * Traverses triplets of all particles in the given cluster. Always uses newton 3 in the AoS data layout.
   * @param cluster The cluster to traverse.
   */
  void traverseClusterTriwise(internal::Cluster<Particle> &cluster) {
    if (_dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = i + 1; j < _clusterSize; j++) {
          for (size_t k = j + 1; k < _clusterSize; k++) {
            if (_useNewton3) {
              _functor->AoSFunctor(cluster[i], cluster[j], cluster[k], true);
            } else {
              _functor->AoSFunctor(cluster[i], cluster[j], cluster[k], false);
              _functor->AoSFunctor(cluster[j], cluster[k], cluster[i], false);
              _functor->AoSFunctor(cluster[k], cluster[i], cluster[j], false);
            }
          }
        }
      }
    } else {
      _functor->SoAFunctorSingle(cluster.getSoAView(), _useNewton3);
    }
  }

  /**
   * Traverses all pairs of particles between two clusters.
   * @param cluster The first cluster.
   * @param neighborCluster The second cluster.
   */
  void traverseClusterPair(internal::Cluster<Particle> &cluster, internal::Cluster<Particle> &neighborCluster) {
    if (_dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = 0; j < _clusterSize; j++) {
          _functor->AoSFunctor(cluster[i], neighborCluster[j], _useNewton3);
        }
      }
    } else {
      _functor->SoAFunctorPair(cluster.getSoAView(), neighborCluster.getSoAView(), _useNewton3);
    }
  }

  /**
   * Traverses all triplets of particles between three clusters.
   * @param cluster The first cluster.
   * @param neighborCluster1 The first neighbor cluster.
   * @param neighborCluster2 The second neighbor cluster.
   */
  void traverseClusterTriplet(internal::Cluster<Particle> &cluster, internal::Cluster<Particle> &neighborCluster1, internal::Cluster<Particle> &neighborCluster2) {
    if (_dataLayout == DataLayoutOption::aos){
      // 3 remaining options for particle triplets:
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = 0; j < _clusterSize; j++) {
          for (size_t k = 0; k < _clusterSize; k++) {
            // 2 particles in cluster + 1 particle in neighbourCluster1 OR neighbourCluster2
            if ((!_useNewton3 && i != j) || i < j) {
              _functor->AoSFunctor(cluster[i], cluster[j], neighborCluster1[k], _useNewton3);
              _functor->AoSFunctor(cluster[i], cluster[j], neighborCluster2[k], _useNewton3);
            }
            // 1 particle in cluster + 2 particles both in neighbourCluster1 OR both in neighbourCluster2
            if (j < k){
              _functor->AoSFunctor(cluster[i], neighborCluster1[j], neighborCluster1[k], _useNewton3);
              _functor->AoSFunctor(cluster[i], neighborCluster2[j], neighborCluster2[k], _useNewton3);
            }
            // 1 particle in cluster + 1 particle in neighbourCluster1 AND neighbourCluster2 each
            _functor->AoSFunctor(cluster[i], neighborCluster1[j], neighborCluster2[k], _useNewton3);
          }
        }
      }
    }
    else {
      _functor->SoAFunctorTriple(cluster.getSoAView(), neighborCluster1.getSoAView(), neighborCluster2.getSoAView(), _useNewton3);
    }
  }

 private:
  Functor *_functor;
  size_t _clusterSize;
  DataLayoutOption _dataLayout;
  bool _useNewton3;
};

}  // namespace autopas::internal
