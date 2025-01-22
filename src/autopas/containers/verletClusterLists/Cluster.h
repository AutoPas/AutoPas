/**
 * @file Cluster.h
 * @author humig
 * @date 27.07.19
 */

#pragma once

#include <algorithm>
#include <iterator>
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
   * @param clusterSize Number of particles in the cluster.
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
   * Indicates if the cluster contains any non-dummy particles.
   * @return True if only dummy particles are present in the cluster.
   */
  bool empty() const {
    const auto firstNonDummy = std::find_if(_firstParticle, _firstParticle + _clusterSize,
                                            [](const auto &particle) { return not particle.isDummy(); });
    return firstNonDummy > &_firstParticle[_clusterSize - 1];
  }

  /**
   * Get Minimum and Maximum of the particles in z-direction.
   * @note This assumes that the particles are sorted along the z-direction!
   * @return Tuple of minimum and maximum in z-direction and bool indicating whether this cluster contains at least one
   * proper particle.
   */
  [[nodiscard]] std::tuple<double, double, bool> getZMinMax() const {
    // Find first particle which is not a dummy and get its z-value.
    auto *begin = &operator[](0);
    auto *end = &operator[](_clusterSize);
    auto pfirst = std::find_if(begin, end, [](const auto &particle) { return not particle.isDummy(); });
    double min = pfirst != end ? pfirst->getR()[2] : std::numeric_limits<double>::max();

    // Find last particle which is not a dummy and get its z-value.
    auto rbegin = std::make_reverse_iterator(end);
    auto rend = std::make_reverse_iterator(begin);
    auto plast = std::find_if(rbegin, rend, [](const auto &particle) { return not particle.isDummy(); });
    double max = plast != rend ? plast->getR()[2] : std::numeric_limits<double>::min();

    return {min, max, pfirst != end};
  }

  /**
   * Returns the SoAView for this cluster.
   * @return the SoAView for this cluster.
   */
  auto getSoAView() { return _soaView; }

  /**
   * Set the SoAView for this cluster.
   * @param view the new SoAView for this cluster.
   */
  void setSoAView(const SoAView<typename Particle::SoAArraysType> &view) { _soaView = view; }

  /**
   * Set the internal neighbor list pointer to an allocated, but not necessarily complete, existing list.
   * @param neighborList Allocated neighbor list.
   */
  void setNeighborList(std::vector<Cluster<Particle> *> *neighborList) { _neighborClusters = neighborList; }

  /**
   * Set the internal pair neighbor list pointer to an allocated, but not necessarily complete, existing list.
   * @param pairNeighborList Allocated pair neighbor list.
   */
  void setPairNeighborList(std::vector<std::pair<Cluster<Particle> *, Cluster<Particle> *>> *pairNeighborList) { _pairNeighborClusters = pairNeighborList; }

  /**
   * Returns the reference to the neighbor list for this cluster.
   * @return reference to the neighbor list.
   */
  std::vector<Cluster<Particle> *> *getNeighbors() { return _neighborClusters; }

  /**
   * Returns the reference to the pair neighbor list for this cluster.
   * @return reference to the pair neighbor list.
   */
  std::vector<std::pair<Cluster *, Cluster *>> *getNeighborPairs() { return _pairNeighborClusters; }

  /**
   * Adds the given cluster to the neighbor list of this cluster.
   * @param neighbor The cluster to add as neighbor.
   */
  void addNeighbor(Cluster<Particle> &neighbor) { _neighborClusters->push_back(&neighbor); }

  /**
   * Adds the two neighbor clusters to the neighbor list of this cluster as a pair.
   * @param neighbor1 First cluster of the neighbor pair
   * @param neighbor2 Second cluster of the neighbor pair
   */
  void addNeighborPair(Cluster<Particle> *neighbor1, Cluster<Particle> *neighbor2) {
    const auto neighborPair = std::make_pair(neighbor1, neighbor2);
    _pairNeighborClusters->push_back(neighborPair);
  }

  /**
   * Remove all neighbors.
   */
  void clearNeighbors() {
    if (_neighborClusters) {
      _neighborClusters->clear();
    }
  }

  /**
   * Remove all pair neighbors.
   */
  void clearPairNeighbors() {
    if (_pairNeighborClusters) {
      _pairNeighborClusters->clear();
    }
  }

  /**
   *
   * @param firstParticle
   */
  void reset(Particle *firstParticle) { _firstParticle = firstParticle; }

  /**
   * Get the bounding box of this cluster
   * @return tuple<lowerCorner, upperCorner>
   */
  [[nodiscard]] std::tuple<std::array<double, 3>, std::array<double, 3>> getBoundingBox() const {
    auto lowerCorner = _firstParticle->getR();
    auto upperCorner = _firstParticle[_clusterSize - 1].getR();

    for (size_t i = 0; i < _clusterSize; ++i) {
      const auto &pos = _firstParticle[i].getR();
      // no need to check z direction, this is already correct by initialization thanks to the particles being sorted in
      // z dimension.
      for (size_t dim = 0; dim < 2; ++dim) {
        lowerCorner[dim] = std::min(lowerCorner[dim], pos[dim]);
        upperCorner[dim] = std::max(upperCorner[dim], pos[dim]);
      }
    }
    return {lowerCorner, upperCorner};
  }

 private:
  /**
   * The number of particles in a full cluster.
   */
  size_t _clusterSize;

  /**
   * A pointer to the first particle of the cluster.
   */
  Particle *_firstParticle = nullptr;
  /**
   * The SoAView for this cluster.
   */
  SoAView<typename Particle::SoAArraysType> _soaView;
  /**
   * The list of neighbor clusters of this cluster.
   */
  std::vector<Cluster *> *_neighborClusters = nullptr;
  /**
   * The list of pair neighbor clusters of this cluster. Only needed for specific 3b-traversals.
   */
  std::vector<std::pair<Cluster *, Cluster *>> *_pairNeighborClusters = nullptr;
};

}  // namespace autopas::internal
