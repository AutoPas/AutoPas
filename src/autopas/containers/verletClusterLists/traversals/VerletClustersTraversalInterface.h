/**
 * @file VerletClustersTraversalInterface.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "../VerletClusterMaths.h"

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
   * @param cellsPerDim the cells per dimension of the container.
   * @param clusterSize the cluster size.
   * @param clusters the clusters of the container.
   * @param neighborLists the neighbor lists of the container.
   */
  virtual void traverseParticlePairs(std::array<VerletClusterMaths::index_t, 3> cellsPerDim, int clusterSize,
                                     std::vector<FullParticleCell<Particle>> &clusters,
                                     std::vector<std::vector<std::vector<Particle *>>> &neighborLists) = 0;
};
}  // namespace autopas
