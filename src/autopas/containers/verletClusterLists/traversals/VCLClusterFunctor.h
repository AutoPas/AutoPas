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
          Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
          // 2 options involving 2 clusters:
          // - 1 particle in cluster + 2 particles in neighbor cluster
          // - 2 particles in cluster + 1 particle in neighbor cluster
          traverseClusterPairTriwise(cluster, neighbor1);
          for (auto neighborClusterIter2 = neighborClusterIter1 + 1; neighborClusterIter2 != neighborClustersEnd; ++neighborClusterIter2) {
            Cluster<Particle> &neighbor2 = **(neighborClusterIter2);
            // 1 option involving all 3 clusters:
            // - one particle in cluster, neighbor cluster 1 and neighbor cluster 2 each
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
    //TODO: switch depending on data layout
    if (not isHaloCluster) {
      traverseClusterTriwise(cluster);
    }
    if (*(cluster.getNeighbors()->data())) {
      auto intersectingNeighbors = std::vector<internal::Cluster<Particle>*>();

      auto &neighborClusters1 = *(cluster.getNeighbors());
      std::sort(neighborClusters1.begin(), neighborClusters1.end());
      auto neighborClusters1End = neighborClusters1.end();

      for (auto neighborClusterIter1 = neighborClusters1.begin(); neighborClusterIter1 != neighborClusters1End; ++neighborClusterIter1) {
        Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
        //take care of all neighbor pairs first
        traverseClusterPairTriwise(cluster, neighbor1);
        auto &neighborClusters2 = *(neighbor1.getNeighbors());
        std::sort(neighborClusters2.begin(), neighborClusters2.end());
        //reserve space in buffer
        std::size_t maxIntersectionSize = std::min(neighborClusters1.size(), neighborClusters2.size());
        intersectingNeighbors.reserve(maxIntersectionSize);
        //intersect the neighbors of cluster and neighbor1
        auto intersectingNeighborsIter = neighborClusterIter1;
        auto intersectingNeighborsEnd = std::set_intersection(++intersectingNeighborsIter, neighborClusters1End,
                                                             neighborClusters2.begin(), neighborClusters2.end(),
                                                             std::back_inserter(intersectingNeighbors));
        //iterate over intersected list
        for (auto neighborClusterIter2 = intersectingNeighbors.begin(); neighborClusterIter2 != intersectingNeighbors.end(); ++neighborClusterIter2) {
          Cluster<Particle> &neighbor2 = **(neighborClusterIter2);
          traverseClusterTriplet(cluster, neighbor1, neighbor2);
        }
        //clear buffer
        intersectingNeighbors.clear();
      }
    }
  }

  /**
   * Invokes the calculations of all interactions within the cluster and, if they exist, it's neighbors using pair lists
   * @param cluster
   * @param isHaloCluster
   */
  void processClusterPairListInteration(internal::Cluster<Particle> &cluster, bool isHaloCluster) {
    //TODO: switch depending on data layout
    if (not isHaloCluster) {
      traverseClusterTriwise(cluster);
    }
    if (*(cluster.getNeighbors()->data())) {
      auto &neighborClusters = *(cluster.getNeighbors());
      auto neighborClustersEnd = neighborClusters.end();
      for (auto neighborClusterIter1 = neighborClusters.begin(); neighborClusterIter1 != neighborClustersEnd; ++neighborClusterIter1) {
        Cluster<Particle> &neighbor1 = **(neighborClusterIter1);
        // 2 options involving 2 clusters:
        // - 1 particle in cluster + 2 particles in neighbor cluster
        // - 2 particles in cluster + 1 particle in neighbor cluster
        traverseClusterPairTriwise(cluster, neighbor1);
      }
    }
    if (cluster.getNeighborPairs()->data()) {
      //traverse pair neighbor list for interactions involving 3 clusters
      auto &neighborClusterPairs = *(cluster.getNeighborPairs());
      auto neighborClusterPairsEnd = neighborClusterPairs.end();
      for (auto neighborClusterPairIter = neighborClusterPairs.begin();
           neighborClusterPairIter != neighborClusterPairsEnd; ++neighborClusterPairIter) {
        Cluster<Particle> &neighbor1 = *((*neighborClusterPairIter).first);
        Cluster<Particle> &neighbor2 = *((*neighborClusterPairIter).second);
        // 1 option involving all 3 clusters:
        // - one particle in cluster, neighbor cluster 1 and neighbor cluster 2 each
        traverseClusterTriplet(cluster, neighbor1, neighbor2);
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
   * Traverses all triplets of particles between two clusters.
   * @param cluster The first cluster.
   * @param neighborCluster The second cluster.
   */
  void traverseClusterPairTriwise(internal::Cluster<Particle> &cluster, internal::Cluster<Particle> &neighborCluster) {
    if (_dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = 0; j < _clusterSize; j++) {
          for (size_t k = j + 1; k < _clusterSize; k++) {
            // 2 particles in cluster + 1 particle in neighbourCluster
            _functor->AoSFunctor(cluster[j], cluster[k], neighborCluster[i], _useNewton3);
            _functor->AoSFunctor(cluster[k], cluster[j], neighborCluster[i], _useNewton3);
            // 1 particle in cluster + 2 particles in neighbourCluster
            _functor->AoSFunctor(cluster[i], neighborCluster[j], neighborCluster[k], _useNewton3);
          }
        }
      }
    } else {
      _functor->SoAFunctorPair(cluster.getSoAView(), neighborCluster.getSoAView(), _useNewton3);
    }
  }

  /**
   * Traverses all triplets of particles between three clusters. All clusters must be part of triplet
   * (Does not work with N3 yet)
   * @param cluster The first cluster.
   * @param neighborCluster1 The first neighbor cluster.
   * @param neighborCluster2 The second neighbor cluster.
   */
  void traverseClusterTriplet(internal::Cluster<Particle> &cluster, internal::Cluster<Particle> &neighborCluster1, internal::Cluster<Particle> &neighborCluster2) {
    if (_dataLayout == DataLayoutOption::aos){
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = 0; j < _clusterSize; j++) {
          for (size_t k = 0; k < _clusterSize; k++) {
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
