/**
 * @file Cluster.h
 * @author humig
 * @date 27.07.19
 */

#pragma once

#include <vector>

#include "autopas/utils/SoAView.h"

namespace autopas::internal {

/**
 * This class represents a cluster in the VerletClusterLists container.
 *
 * It contains a pointer to the particles for AoS, a SoAView for SoA, and the neighbor list for this cluster.
 *
 * @tparam Particle The type of the particles this cluster consists of.
 * @tparam clusterSize The number of particles in the cluster.
 */
template <class Particle>
class Cluster {
 public:
  /**
   * Constructs a cluster starting from firstParticle and going on for clusterSize particles.
   *
   * Caller is responsible that there are enough particles after this particle in memory.
   *
   * @param firstParticle A pointer to the first particle of the cluster.
   */
  explicit Cluster(Particle *firstParticle, size_t clusterSize)
      : _firstParticle(firstParticle), _clusterSize(clusterSize) {}

  /**
   * Returns the particle at position index in the cluster.
   *
   * No index checking is performed!
   *
   * @param index The index of the particle to return.
   * @return the particle at position index in the cluster.
   */
  Particle &operator[](size_t index) { return *(_firstParticle + index); }

  /**
   * @copydoc operator[](size_t)
   */
  const Particle &operator[](size_t index) const { return *(_firstParticle + index); }

  /**
   * Returns the SoAView for this cluster.
   * @return the SoAView for this cluster.
   */
  auto getSoAView() { return _soaView; }

  /**
   * Set the SoAView for this cluster.
   * @param view the new SoAView for this cluster.
   */
  void setSoAView(SoAView<typename Particle::SoAArraysType> view) { _soaView = view; }

  /**
   * Returns the neighbor list for this cluster.
   * @return the neighbor list for this cluster.
   */
  const auto &getNeighbors() const { return _neighborClusters; }

  /**
   * Adds the given cluster to the neighbor list of this cluster.
   * @param neighbor The cluster to add as neighbor.
   */
  void addNeighbor(Cluster<Particle> &neighbor) { _neighborClusters.push_back(&neighbor); }

  /**
   * Remove all neighbors.
   */
  void clearNeighbors() { _neighborClusters.clear(); }

 private:
  /**
   * The number of particles in a full cluster.
   */
  size_t _clusterSize;

  /**
   * A pointer to the first particle of the cluster.
   */
  Particle *_firstParticle;
  /**
   * The SoAView for this cluster.
   */
  SoAView<typename Particle::SoAArraysType> _soaView;
  /**
   * The list of neighbor clusters of this cluster.
   */
  std::vector<Cluster *> _neighborClusters;
};

}  // namespace autopas::internal
