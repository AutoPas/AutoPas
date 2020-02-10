/**
 * @file VerletClustersTraversalInterface.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

namespace autopas {

template <class Particle>
class VerletClusterLists;

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
   * Sets the cluster list for the traversal to iterate over.
   * @param verletClusterLists the cluster list to iterate over.
   */
  virtual void setClusterLists(VerletClusterLists<Particle> &verletClusterLists) {
    _verletClusterLists = &verletClusterLists;
  }

  /**
   * Returns whether this traversal needs the static cluster thread partiton of the cluster list.
   * @returns whether this traversal needs the static cluster thread partiton of the cluster list.
   */
  virtual bool needsStaticClusterThreadPartition() { return false; };

 protected:
  /**
   * The cluster list to iterate over.
   */
  VerletClusterLists<Particle> *_verletClusterLists;
};
}  // namespace autopas
