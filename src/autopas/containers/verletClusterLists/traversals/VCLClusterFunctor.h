/**
 * @file VCLClusterFunctor.h
 * @author humig
 * @date 12.08.19
 */

#pragma once

#include "autopas/containers/verletClusterLists/Cluster.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::internal {
/**
 * Provides methods to traverse a single cluster and a pair of clusters.
 *
 * @tparam Particle_T The type of particle the clusters contain.
 * @tparam PairwiseFunctor The type of the functor the VCLClusterFunctor should use.
 */
template <class Particle_T, class PairwiseFunctor>
class VCLClusterFunctor {
 public:
  /**
   * Constructs a VCLClusterFunctor that uses the given functor internally.
   * @param functor The functor to use internally.
   * @param clusterSize Number of particles per cluster.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit VCLClusterFunctor(PairwiseFunctor *functor, size_t clusterSize, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(functor), _clusterSize(clusterSize), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Invokes the calculations of all interactions within the cluster and, if they exist, it's neighbors.
   * @param cluster
   * @param isHaloCluster
   */
  void processCluster(internal::Cluster<Particle_T> &cluster, bool isHaloCluster) {
    if (not isHaloCluster) {
      traverseCluster(cluster);
    }

    // only iterate neighbors if the neighbor list contains more than just nullptr
    if (*(cluster.getNeighbors()->data())) {
      for (auto *neighborClusterPtr : *(cluster.getNeighbors())) {
        traverseClusterPair(cluster, *neighborClusterPtr);
      }
    }
  }

 private:
  /**
   * Traverses pairs of all particles in the given cluster. Always uses newton 3 in the AoS data layout.
   * @param cluster The cluster to traverse.
   */
  void traverseCluster(internal::Cluster<Particle_T> &cluster) {
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
   * Traverses all pairs of particles between two clusters.
   * @param cluster The first cluster.
   * @param neighborCluster The second cluster.
   */
  void traverseClusterPair(internal::Cluster<Particle_T> &cluster, internal::Cluster<Particle_T> &neighborCluster) {
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

 private:
  PairwiseFunctor *_functor;
  size_t _clusterSize;
  DataLayoutOption _dataLayout;
  bool _useNewton3;
};

}  // namespace autopas::internal
