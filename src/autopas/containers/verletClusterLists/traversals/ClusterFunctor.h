/**
 * @file ClusterFunctor.h
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
 * @tparam PairwiseFunctor The type of the functor the ClusterFunctor should use.
 * @tparam dataLayout The data layout to use.
 * @tparam useNewton3 If newton 3 should be used or not.
 * @tparam clusterSize The number of particles in each cluster.
 */
template <class Particle, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3,
          size_t clusterSize>
class ClusterFunctor {
 public:
  /**
   * Constructs a ClusterFunctor that uses the given functor internally.
   * @param functor The functor to use internally.
   */
  explicit ClusterFunctor(PairwiseFunctor *functor) : _functor(functor) {}

  /**
   * Traverses pairs of all particles in the given cluster. Always uses newton 3 in the AoS data layout.
   * @param cluster The cluster to traverse.
   */
  void traverseCluster(internal::Cluster<Particle, clusterSize> &cluster) {
    if constexpr (dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < clusterSize; i++) {
        // Always use newton 3 for interactions within one cluster.
        for (size_t j = i + 1; j < clusterSize; j++) {
          _functor->AoSFunctor(cluster[i], cluster[j], true);
        }
      }
    } else {
      _functor->SoAFunctor(cluster.getSoAView(), useNewton3);
    }
  }

  /**
   * Traverses all pairs of particles between two clusters.
   * @param cluster The first cluster.
   * @param neighborCluster The second cluster.
   */
  void traverseClusterPair(internal::Cluster<Particle, clusterSize> &cluster,
                           internal::Cluster<Particle, clusterSize> &neighborCluster) {
    if constexpr (dataLayout == DataLayoutOption::aos) {
      for (size_t i = 0; i < clusterSize; i++) {
        for (size_t j = 0; j < clusterSize; j++) {
          _functor->AoSFunctor(cluster[i], neighborCluster[j], useNewton3);
        }
      }
    } else {
      _functor->SoAFunctor(cluster.getSoAView(), neighborCluster.getSoAView(), useNewton3);
    }
  }

 private:
  PairwiseFunctor *_functor;
};

}  // namespace autopas::internal
