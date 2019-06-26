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
   * Sets the information that traversals over the verlet cluster lists container can use.
   * @param verletClusterLists The container to traverse.
   */
  virtual void setTraversalInfo(VerletClusterLists<Particle> *verletClusterLists) {
    _verletClusterLists = verletClusterLists;
  }

  /**
   * Iterates over all particle pairs.
   */
  virtual void traverseParticlePairs() = 0;

  /**
   * Initializes the traversal of the container.
   */
  virtual void initClusterTraversal() = 0;

  /**
   * Finalizes the traversal over the container.
   */
  virtual void endClusterTraversal() = 0;

 protected:
  /// The container to traverse.
  VerletClusterLists<Particle> *_verletClusterLists;
};
}  // namespace autopas
