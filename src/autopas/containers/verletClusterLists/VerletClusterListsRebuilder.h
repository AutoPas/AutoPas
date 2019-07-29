/**
 * @file Cluster.h
 * @author humig
 * @date 29.07.19
 */

#pragma once

#include "VerletClusterLists.h"
#include "autopas/utils/inBox.h"

namespace autopas {

template <class Particle>
class VerletClusterLists;

namespace internal {

template <class Particle>
class VerletClusterListsRebuilder {
 private:
  static constexpr size_t clusterSize = VerletClusterLists<Particle>::clusterSize;

  std::vector<Particle> &_particlesToAdd;
  std::vector<ClusterTower<Particle, clusterSize>> &_towers;
  double _towerSideLength;
  int _interactionLengthInTowers;
  double _towerSideLengthReciprocal;
  std::array<size_t, 3> _towersPerDim;
  bool _neighborListIsNewton3;
  double _interactionLength;
  double _interactionLengthSqr;
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;

 public:
  VerletClusterListsRebuilder(VerletClusterLists<Particle> &clusterList, std::vector<Particle> &particlesToAdd,
                              bool newton3)
      : _particlesToAdd(particlesToAdd),
        _towers(clusterList.getTowers()),
        _towerSideLength(clusterList.getTowerSideLength()),
        _interactionLengthInTowers(clusterList.getInteractionLengthInTowers()),
        _towerSideLengthReciprocal(1 / _towerSideLength),
        _towersPerDim(clusterList.getTowersPerDimension()),
        _neighborListIsNewton3(newton3),
        _interactionLength(clusterList.getInteractionLength()),
        _interactionLengthSqr(_interactionLength * _interactionLength),
        _boxMin(clusterList.getBoxMin()),
        _boxMax(clusterList.getBoxMax()) {}

  struct BuildingReturnValues {
    double _towerSideLength;
    int _interactionLengthInTowers;
    std::array<size_t, 3> _towersPerDim;
    size_t _numClusters;
    bool _neighborListIsNewton3;
  };
  /**
   * Recalculate towers and clusters, build neighbor lists and pad clusters.
   *
   * @returns tuple consisting of
   */
  BuildingReturnValues rebuild() {
    std::vector<Particle> invalidParticles = collectAllParticlesFromTowers();
    invalidParticles.insert(invalidParticles.end(), _particlesToAdd.begin(), _particlesToAdd.end());
    _particlesToAdd.clear();

    auto boxSize = ArrayMath::sub(_boxMax, _boxMin);

    _towerSideLength = estimateOptimalGridSideLength(invalidParticles.size(), boxSize);
    _interactionLengthInTowers = static_cast<int>(std::ceil(_interactionLength / _towerSideLength));
    _towerSideLengthReciprocal = 1 / _towerSideLength;

    _towersPerDim = calculateTowersPerDim(boxSize);
    // _towersPerDim[2] is always 1
    size_t numTowers = _towersPerDim[0] * _towersPerDim[1];

    // resize to number of towers
    _towers.resize(numTowers);

    sortParticlesIntoTowers(invalidParticles);

    size_t numClusters = 0;
    for (auto &tower : _towers) {
      numClusters += tower.generateClusters();
    }

    updateNeighborLists();

    // Make sure that each cluster contains a multiple of clusterSize particles, and that dummies don't interact.
    double dummyParticleDistance = _interactionLength * 2;
    double startDummiesX = 1000 * _boxMax[0];
    for (size_t index = 0; index < _towers.size(); index++) {
      _towers[index].fillUpWithDummyParticles(startDummiesX + index * dummyParticleDistance, dummyParticleDistance);
    }

    return {_towerSideLength, _interactionLengthInTowers, _towersPerDim, numClusters, _neighborListIsNewton3};
  }

 protected:
  /**
   * Takes all particles from all towers and returns them. Towers are cleared afterwards.
   * @return All particles in the container.
   */
  std::vector<Particle> collectAllParticlesFromTowers() {
    std::vector<Particle> invalidParticles;
    for (auto &tower : _towers) {
      invalidParticles.reserve(invalidParticles.size() + tower.getNumActualParticles());
      for (auto it = tower.begin(); it.isValid(); ++it) {
        invalidParticles.push_back(*it);
      }
      tower.clear();
    }
    return invalidParticles;
  }

  /**
   * Estimates the optimal grid side length.
   * @param numParticles The number of particles in the container.
   * @param boxSize The size of the domain.
   * @return an estimated optimal grid side length.
   */
  [[nodiscard]] virtual double estimateOptimalGridSideLength(size_t numParticles, std::array<double, 3> boxSize) const {
    double volume = boxSize[0] * boxSize[1] * boxSize[2];
    if (numParticles > 0) {
      // estimate particle density
      double density = numParticles / volume;

      return std::cbrt(clusterSize / density);
    } else {
      return std::max(boxSize[0], boxSize[1]);
    }
  }

  /**
   * Calculates the cells per dimension in the container using the _towerSideLengthReciprocal.
   * @param boxSize the size of the domain.
   * @return the cells per dimension in the container.
   */
  [[nodiscard]] std::array<size_t, 3> calculateTowersPerDim(std::array<double, 3> boxSize) const {
    std::array<size_t, 3> towersPerDim{};
    for (int d = 0; d < 2; d++) {
      towersPerDim[d] = static_cast<size_t>(std::ceil(boxSize[d] * _towerSideLengthReciprocal));
      // at least one cell
      towersPerDim[d] = std::max(towersPerDim[d], 1ul);
    }
    towersPerDim[2] = 1ul;
    return towersPerDim;
  }

  /**
   * Sorts all passed particles in the appropriate clusters.
   * @param particles The particles to sort in the clusters.
   */
  void sortParticlesIntoTowers(const std::vector<Particle> &particles) {
    const auto numParticles = particles.size();
#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic, 1000) default(none) shared(particles, numParticles)
#endif
    for (size_t index = 0; index < numParticles; index++) {
      const Particle &particle = particles[index];
      if (utils::inBox(particle.getR(), _boxMin, _boxMax)) {
        auto &tower = getTowerForParticleLocation(particle.getR());
        tower.addParticle(particle);
      }
    }
  }

  /**
   * Updates the neighbor lists.
   */
  void updateNeighborLists() {
    const int maxTowerIndexX = _towersPerDim[0] - 1;
    const int maxTowerIndexY = _towersPerDim[1] - 1;

    // for all towers
#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2) default(none) shared(maxTowerIndexX, maxTowerIndexY)
#endif
    for (int towerIndexY = 0; towerIndexY <= maxTowerIndexY; towerIndexY++) {
      for (int towerIndexX = 0; towerIndexX <= maxTowerIndexX; towerIndexX++) {
        const int minX = std::max(towerIndexX - _interactionLengthInTowers, 0);
        const int minY = std::max(towerIndexY - _interactionLengthInTowers, 0);
        const int maxX = std::min(towerIndexX + _interactionLengthInTowers, maxTowerIndexX);
        const int maxY = std::min(towerIndexY + _interactionLengthInTowers, maxTowerIndexY);

        calculateNeighborsForTowerInRange(towerIndexX, towerIndexY, minX, maxX, minY, maxY);
      }
    }
  }

  /**
   * Calculates the neighbors for all clusters of the given tower that are in the given neighbor towers.
   * @param towerIndexX The x-coordinate of the given tower.
   * @param towerIndexY the y-coordinate of the given tower.
   * @param minX The minimum x-coordinate of the given neighbor towers.
   * @param maxX The maximum x-coordinate of the given neighbor towers.
   * @param minY The minimum y-coordinate of the given neighbor towers.
   * @param maxY The maximum y-coordinate of the given neighbor towers.
   */
  void calculateNeighborsForTowerInRange(const int towerIndexX, const int towerIndexY, const int minX, const int maxX,
                                         const int minY, const int maxY) {
    auto &tower = getTowerAtCoordinates(towerIndexX, towerIndexY);
    // for all neighbor towers
    for (int neighborIndexY = minY; neighborIndexY <= maxY; neighborIndexY++) {
      double distBetweenTowersY = std::max(0, std::abs(towerIndexY - neighborIndexY) - 1) * _towerSideLength;

      for (int neighborIndexX = minX; neighborIndexX <= maxX; neighborIndexX++) {
        if (_neighborListIsNewton3 and not shouldTowerContainOtherAsNeighborWithNewton3(
                                           towerIndexX, towerIndexY, neighborIndexX, neighborIndexY)) {
          continue;
        }

        double distBetweenTowersX = std::max(0, std::abs(towerIndexX - neighborIndexX) - 1) * _towerSideLength;

        // calculate distance in xy-plane and skip if already longer than cutoff
        auto distBetweenTowersXYsqr = distBetweenTowersX * distBetweenTowersX + distBetweenTowersY * distBetweenTowersY;
        if (distBetweenTowersXYsqr <= _interactionLengthSqr) {
          auto &neighborTower = getTowerAtCoordinates(neighborIndexX, neighborIndexY);

          calculateNeighborsForTowerPair(tower, neighborTower, distBetweenTowersXYsqr);
        }
      }
    }
  }

  int get1DInteractionCellIndexForTower(const int towerIndexX, const int towerIndexY) {
    const int interactionCellTowerX = towerIndexX / _interactionLengthInTowers;
    const int interactionCellTowerY = towerIndexY / _interactionLengthInTowers;

    const int numInteractionCellsX = static_cast<int>(std::ceil(_towersPerDim[0] / (double)_interactionLengthInTowers));

    return interactionCellTowerX + numInteractionCellsX * interactionCellTowerY;
  }

  bool shouldTowerContainOtherAsNeighborWithNewton3(const int towerIndexX, const int towerIndexY,
                                                    const int neighborIndexX, const int neighborIndexY) {
    auto interactionCellTowerIndex1D = get1DInteractionCellIndexForTower(towerIndexX, towerIndexY);
    auto interactionCellNeighborIndex1D = get1DInteractionCellIndexForTower(neighborIndexX, neighborIndexY);

    auto towerIndex1D = towerIndex2DTo1D(towerIndexX, towerIndexY);
    auto neighborIndex1D = towerIndex2DTo1D(neighborIndexX, neighborIndexY);

    return interactionCellNeighborIndex1D > interactionCellTowerIndex1D or
           (interactionCellNeighborIndex1D == interactionCellTowerIndex1D and neighborIndex1D >= towerIndex1D);
  }

  void calculateNeighborsForTowerPair(internal::ClusterTower<Particle, clusterSize> &tower,
                                      internal::ClusterTower<Particle, clusterSize> &neighborTower,
                                      double distBetweenTowersXYsqr) {
    const bool isSameTower = &tower == &neighborTower;
    for (size_t towerIndex = 0; towerIndex < tower.getNumClusters(); towerIndex++) {
      auto startIndexNeighbor = _neighborListIsNewton3 and isSameTower ? towerIndex + 1 : 0;
      auto &towerCluster = tower.getCluster(towerIndex);
      double towerClusterBoxBottom = towerCluster.getParticle(0).getR()[2];
      double towerClusterBoxTop = towerCluster.getParticle(clusterSize - 1).getR()[2];

      for (size_t neighborIndex = startIndexNeighbor; neighborIndex < neighborTower.getNumClusters(); neighborIndex++) {
        const bool isSameCluster = towerIndex == neighborIndex;
        if (not _neighborListIsNewton3 and isSameTower and isSameCluster) {
          continue;
        }
        auto &neighborCluster = neighborTower.getCluster(neighborIndex);
        double neighborClusterBoxBottom = neighborCluster.getParticle(0).getR()[2];
        double neighborClusterBoxTop = neighborCluster.getParticle(clusterSize - 1).getR()[2];

        double distZ =
            bboxDistance(towerClusterBoxBottom, towerClusterBoxTop, neighborClusterBoxBottom, neighborClusterBoxTop);
        if (distBetweenTowersXYsqr + distZ * distZ <= _interactionLengthSqr) {
          towerCluster.addNeighbor(neighborCluster);
        }
      }
    }
  }

  /**
   * Calculates the distance of two bounding boxes in one dimension.
   * @param min1 minimum coordinate of first bbox in tested dimension
   * @param max1 maximum coordinate of first bbox in tested dimension
   * @param min2 minimum coordinate of second bbox in tested dimension
   * @param max2 maximum coordinate of second bbox in tested dimension
   * @return distance
   */
  [[nodiscard]] double bboxDistance(const double min1, const double max1, const double min2, const double max2) const {
    if (max1 < min2) {
      return min2 - max1;
    } else if (min1 > max2) {
      return min1 - max2;
    } else {
      return 0;
    }
  }

  auto &getTowerForParticleLocation(std::array<double, 3> location) {
    std::array<size_t, 2> cellIndex{};

    for (int dim = 0; dim < 2; dim++) {
      const long int value =
          (static_cast<long int>(floor((location[dim] - _boxMin[dim]) * _towerSideLengthReciprocal))) + 1l;
      const size_t nonNegativeValue = static_cast<size_t>(std::max(value, 0l));
      const size_t nonLargerValue = std::min(nonNegativeValue, _towersPerDim[dim] - 1);
      cellIndex[dim] = nonLargerValue;
      // @todo this is a sanity check to prevent doubling of particles, but could be done better! e.g. by border and
      // flag manager
      if (location[dim] >= _boxMax[dim]) {
        cellIndex[dim] = _towersPerDim[dim] - 1;
      } else if (location[dim] < _boxMin[dim]) {
        cellIndex[dim] = 0;
      }
    }

    return getTowerAtCoordinates(cellIndex[0], cellIndex[1]);
  }

  size_t towerIndex2DTo1D(const size_t x, const size_t y) {
    return VerletClusterLists<Particle>::towerIndex2DTo1D(x, y, _towersPerDim);
  }

  auto &getTowerAtCoordinates(const size_t x, const size_t y) { return _towers[towerIndex2DTo1D(x, y)]; }
};

}  //  namespace internal

}  //  namespace autopas