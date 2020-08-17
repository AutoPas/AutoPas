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
 * @tparam Particle The type of particle the clusters contain.
 * @tparam PairwiseFunctor The type of the functor the VCLClusterFunctor should use.
 * @tparam dataLayout The data layout to use.
 * @tparam useNewton3 If newton 3 should be used or not.
 * @tparam clusterSize The number of particles in each cluster.
 */
template <class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VCLClusterFunctor {
 public:
  /**
   * Constructs a VCLClusterFunctor that uses the given functor internally.
   * @param functor The functor to use internally.
   * @param clusterSize Number of particles per cluster.
   */
  explicit VCLClusterFunctor(PairwiseFunctor *functor, size_t clusterSize)
      : _functor(functor), _clusterSize(clusterSize) {}

  /**
   * Traverses pairs of all particles in the given cluster. Always uses newton 3 in the AoS data layout.
   * @param cluster The cluster to traverse.
   */
  void traverseCluster(internal::Cluster<Particle> &cluster) {
    if constexpr (dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = i + 1; j < _clusterSize; j++) {
          // this if else branch is needed because of https://github.com/AutoPas/AutoPas/issues/426
          if constexpr (useNewton3) {
            _functor->AoSFunctor(cluster[i], cluster[j], true);
          } else {
            _functor->AoSFunctor(cluster[i], cluster[j], false);
            _functor->AoSFunctor(cluster[j], cluster[i], false);
          }
        }
      }
    } else {
      _functor->SoAFunctorSingle(cluster.getSoAView(), useNewton3);
    }
  }

  /**
   * Traverses all pairs of particles between two clusters.
   * @param cluster The first cluster.
   * @param neighborCluster The second cluster.
   */
  void traverseClusterPair(internal::Cluster<Particle> &cluster, internal::Cluster<Particle> &neighborCluster) {
    if constexpr (dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < _clusterSize; i++) {
        for (size_t j = 0; j < _clusterSize; j++) {
          _functor->AoSFunctor(cluster[i], neighborCluster[j], useNewton3);
        }
      }
    } else {
      _functor->SoAFunctorPair(cluster.getSoAView(), neighborCluster.getSoAView(), useNewton3);
    }
  }

 private:
  PairwiseFunctor *_functor;
  size_t _clusterSize;
};

}  // namespace autopas::internal
