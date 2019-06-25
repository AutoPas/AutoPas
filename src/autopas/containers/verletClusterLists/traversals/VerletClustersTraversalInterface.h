/**
 * @file VerletClustersTraversalInterface.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/VerletClusterMaths.h"

namespace autopas {

/**
 * Interface for traversals of the VerletClusterLists container.
 */
template <class Particle>
class VerletClustersTraversalInterface {
 public:
  /**
   * virtual default destructor.
   */
  virtual ~VerletClustersTraversalInterface() = default;

  /**
   * Iterates over all particle pairs.
   * @param verletClusterLists The container to traverse.
   */
  virtual void traverseParticlePairs(VerletClusterLists<Particle> &verletClusterLists) = 0;
};
}  // namespace autopas
