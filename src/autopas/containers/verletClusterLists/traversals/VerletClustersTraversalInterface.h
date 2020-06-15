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
   * Sets the towers of the cluster list for the traversal to iterate over.
   * @param towers towers of the cluster list for the traversal to iterate over.
   */
  virtual void setTowers(std::vector<internal::ClusterTower<Particle>> &towers) { _towers = &towers; }
  /**
   * Returns whether this traversal needs the static cluster thread partiton of the cluster list.
   * @returns whether this traversal needs the static cluster thread partiton of the cluster list.
   */
  virtual bool needsStaticClusterThreadPartition() { return false; };

 protected:
  /**
   * The cluster list to iterate over.
   *
   * It provides methods to iterate over the clusters. If more control over the iteration is needed, traversals can also
   * iterate over the towers directly. They are accessible through _towers.
   */
  VerletClusterLists<Particle> *_verletClusterLists;

  /**
   * The towers of the cluster list to iterate over. These directly contain the particles to be modified by the
   * traversal.
   */
  std::vector<internal::ClusterTower<Particle>> *_towers;
};
}  // namespace autopas
