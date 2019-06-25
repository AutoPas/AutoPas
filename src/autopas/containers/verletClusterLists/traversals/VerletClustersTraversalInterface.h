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
   * Sets the information that traversals over the verlet cluster lists container can use.
   * @param cellsPerDim the cells per dimension of the container.
   * @param numClustesr The number of clusters in the container.
   * @param clusterSize the cluster size.
   * @param grids the grids of the container.
   * @param neighborLists the neighbor lists of the container.
   */
  virtual void setTraversalInfo(std::array<VerletClusterMaths::index_t, 3> cellsPerDim,
                                VerletClusterMaths::index_t numClusters, int clusterSize,
                                std::vector<FullParticleCell<Particle>> &grids,
                                std::vector<std::vector<std::vector<Particle *>>> &neighborLists,
                                std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    _cellsPerDim = cellsPerDim;
    _numClusters = numClusters;
    _clusterSize = clusterSize;
    _grids = &grids;
    _neighborLists = &neighborLists;
    _aosToSoaMap = &aosToSoaMap;
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
  /// The cells per dimension of the container.
  std::array<VerletClusterMaths::index_t, 3> _cellsPerDim;
  /// The number of clusters in the container.
  VerletClusterMaths::index_t _numClusters;
  /// The number of particles in one cluster in the container.
  int _clusterSize;
  /// The grids of the container.
  std::vector<FullParticleCell<Particle>> *_grids;
  /// The neighbor lists of the container.
  std::vector<std::vector<std::vector<Particle *>>> *_neighborLists;
  /// The map of the container that maps cluster start pointers to the index of the cluster.
  std::unordered_map<Particle *, VerletClusterMaths::index_t> *_aosToSoaMap;
};
}  // namespace autopas
