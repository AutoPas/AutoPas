/**
 * @file Cluster.h
 * @author humig
 * @date 29.07.19
 */

#pragma once

#include "VerletClusterLists.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/inBox.h"

namespace autopas {

template <class Particle>
class VerletClusterLists;

namespace internal {

/**
 * Helper class for rebuilding the VerletClusterLists container.
 * @tparam Particle The type of the particle the container contains.
 */
template <class Particle>
class VerletClusterListsRebuilder {
 private:
  size_t _clusterSize;

  std::vector<Particle> &_particlesToAdd;
  std::vector<ClusterTower<Particle>> &_towers;
  double _towerSideLength;
  int _interactionLengthInTowers;
  double _towerSideLengthReciprocal;
  std::array<size_t, 2> _towersPerDim;
  double _interactionLength;
  double _interactionLengthSqr;
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  std::array<double, 3> _haloBoxMin;
  std::array<double, 3> _haloBoxMax;

 public:
  /**
   * Constructs the builder from the cluster list.
   *
   * @param clusterList The cluster list to rebuild the neighbor lists for.
   * @param towers The towers from the cluster list to rebuild.
   * @param particlesToAdd New particles to add.
   * @param clusterSize Size of the clusters in particles.
   */
  VerletClusterListsRebuilder(const VerletClusterLists<Particle> &clusterList,
                              std::vector<ClusterTower<Particle>> &towers, std::vector<Particle> &particlesToAdd,
                              size_t clusterSize)
      : _clusterSize(clusterSize),
        _particlesToAdd(particlesToAdd),
        _towers(towers),
        _towerSideLength(clusterList.getTowerSideLength()),
        _interactionLengthInTowers(clusterList.getNumTowersPerInteractionLength()),
        _towerSideLengthReciprocal(clusterList.getTowerSideLengthReciprocal()),
        _towersPerDim(clusterList.getTowersPerDimension()),
        _interactionLength(clusterList.getInteractionLength()),
        _interactionLengthSqr(_interactionLength * _interactionLength),
        _boxMin(clusterList.getBoxMin()),
        _boxMax(clusterList.getBoxMax()),
        _haloBoxMin(clusterList.getHaloBoxMin()),
        _haloBoxMax(clusterList.getHaloBoxMax()) {}

  /**
   * Rebuilds the towers, clusters, and neighbor lists.
   *
   * @return new values for VerletClusterLists member variables. They are returned as tuple consisting of:
   * {
   *  double:                The new side length of each tower in xy-direction,
   *  int:                   The interaction length in towers using the new tower side length,
   *  std::array<size_t, 2>: The number of towers in each dimension using the new tower side length,
   *  size_t:                The new number of clusters in the container,
   * }
   */
  auto rebuildTowersAndClusters() {
    auto invalidParticles = collectAllParticlesFromTowers();
    invalidParticles.push_back(std::move(_particlesToAdd));
    _particlesToAdd.clear();

    // count particles by accumulating tower sizes
    size_t numParticles = std::accumulate(std::begin(invalidParticles), std::end(invalidParticles), 0,
                                          [](auto acc, auto &v) { return acc + v.size(); });

    auto boxSizeWithHalo = utils::ArrayMath::sub(_haloBoxMax, _haloBoxMin);

    _towerSideLength = estimateOptimalGridSideLength(numParticles, boxSizeWithHalo);
    _towerSideLengthReciprocal = 1 / _towerSideLength;
    _interactionLengthInTowers = static_cast<int>(std::ceil(_interactionLength * _towerSideLengthReciprocal));

    _towersPerDim = calculateTowersPerDim(boxSizeWithHalo);
    size_t numTowers = _towersPerDim[0] * _towersPerDim[1];

    // resize to number of towers. Cannot use resize since towers are not default constructable.
    _towers.clear();
    _towers.reserve(numTowers);

    // create towers and make an estimate for how many particles memory needs to be allocated
    // 2.7 seems high but gave the best performance when testing
    const size_t sizeEstimation = (static_cast<double>(numParticles) / numTowers) * 2.7;
    for (int i = 0; i < numTowers; ++i) {
      _towers.emplace_back(ClusterTower<Particle>(_clusterSize));
      _towers[i].reserve(sizeEstimation);
    }

    sortParticlesIntoTowers(invalidParticles);

    // generate clusters and count them
    size_t numClusters = 0;
    for (auto &tower : _towers) {
      numClusters += tower.generateClusters();
    }

    return std::make_tuple(_towerSideLength, _interactionLengthInTowers, _towersPerDim, numClusters);
  }

  /**
   * Rebuilds the neighbor lists and fills Clusters with dummies as described in
   * ClusterTower::fillUpWithDummyParticles.
   * @param useNewton3 Specifies, whether neighbor lists should use newton3. This changes the way what the lists
   * contain. If an cluster A interacts with cluster B, then this interaction will either show up only once in the
   * interaction lists of the custers (for newton3 == true) or show up in the interaction lists of both (for newton3 ==
   * false)
   */
  void rebuildNeighborListsAndFillClusters(bool useNewton3) {
    clearNeighborListsAndResetDummies();
    updateNeighborLists(useNewton3);

    double dummyParticleDistance = _interactionLength * 2;
    double startDummiesX = 1000 * _haloBoxMax[0];
    for (size_t index = 0; index < _towers.size(); index++) {
      _towers[index].fillUpWithDummyParticles(startDummiesX + index * dummyParticleDistance, dummyParticleDistance);
    }
  }

  /**
   * Removes previously saved neighbors from clusters and sets the positions of the dummy particles to inside of the
   * cluster. The latter reduces the amount of calculated interaction partners.
   */
  void clearNeighborListsAndResetDummies() {
    for (auto &tower : _towers) {
      tower.setDummyParticlesToLastActualParticle();
      for (auto &cluster : tower.getClusters()) {
        cluster.clearNeighbors();
      }
    }
  }

  /**
   * Takes all particles from all towers and returns them. Towers are cleared afterwards.
   * @return All particles in the container sorted in 2D as they were in the towers.
   */
  std::vector<std::vector<Particle>> collectAllParticlesFromTowers() {
    std::vector<std::vector<Particle>> invalidParticles;
    invalidParticles.resize(_towers.size());
    for (size_t towerIndex = 0; towerIndex < _towers.size(); towerIndex++) {
      invalidParticles[towerIndex] = _towers[towerIndex].collectAllActualParticles();
      _towers[towerIndex].clear();
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

      return std::cbrt(_clusterSize / density);
    } else {
      return std::max(boxSize[0], boxSize[1]);
    }
  }

  /**
   * Calculates the cells per dimension in the container using the _towerSideLengthReciprocal.
   * @param boxSize the size of the domain.
   * @return the cells per dimension in the container.
   */
  [[nodiscard]] std::array<size_t, 2> calculateTowersPerDim(std::array<double, 3> boxSize) const {
    std::array<size_t, 2> towersPerDim{};
    for (int dimension = 0; dimension < 2; dimension++) {
      towersPerDim[dimension] = static_cast<size_t>(std::ceil(boxSize[dimension] * _towerSideLengthReciprocal));
      // at least one cell
      towersPerDim[dimension] = std::max(towersPerDim[dimension], 1ul);
    }
    return towersPerDim;
  }

  /**
   * Sorts all passed particles in the appropriate clusters.
   *
   * @note This Function takes a 2D vector because it expects the layout from the old clusters.
   * The information however, is not utilized hence when in doubt all particles can go in one vector.
   *
   * @param particles2D The particles to sort in the towers.
   */
  void sortParticlesIntoTowers(const std::vector<std::vector<Particle>> &particles2D) {
    const auto numVectors = particles2D.size();
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t index = 0; index < numVectors; index++) {
      const std::vector<Particle> &vector = particles2D[index];
      for (const auto &particle : vector) {
        if (utils::inBox(particle.getR(), _haloBoxMin, _haloBoxMax)) {
          auto &tower = getTower(particle.getR());
          tower.addParticle(particle);
        } else {
          AutoPasLog(trace, "Not adding particle to VerletClusterLists container, because it is far outside:\n{}",
                     particle.toString());
        }
      }
    }
  }

  /**
   * Updates the neighbor lists.
   * @param useNewton3 Specifies, whether neighbor lists should use newton3. This changes the way what the lists
   * contain. If an cluster A interacts with cluster B, then this interaction will either show up only once in the
   * interaction lists of the custers (for newton3 == true) or show up in the interaction lists of both (for newton3 ==
   * false)
   */
  void updateNeighborLists(bool useNewton3) {
    const int maxTowerIndexX = _towersPerDim[0] - 1;
    const int maxTowerIndexY = _towersPerDim[1] - 1;

    // for all towers
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (int towerIndexY = 0; towerIndexY <= maxTowerIndexY; towerIndexY++) {
      for (int towerIndexX = 0; towerIndexX <= maxTowerIndexX; towerIndexX++) {
        const int minX = std::max(towerIndexX - _interactionLengthInTowers, 0);
        const int minY = std::max(towerIndexY - _interactionLengthInTowers, 0);
        const int maxX = std::min(towerIndexX + _interactionLengthInTowers, maxTowerIndexX);
        const int maxY = std::min(towerIndexY + _interactionLengthInTowers, maxTowerIndexY);

        iterateNeighborTowers(towerIndexX, towerIndexY, minX, maxX, minY, maxY, useNewton3,
                              [this](auto &towerA, auto &towerB, double distBetweenTowersXYsqr, bool useNewton3) {
                                calculateNeighborsBetweenTowers(towerA, towerB, distBetweenTowersXYsqr, useNewton3);
                              });
      }
    }
  }

  /**
   * For all clusters in a tower, given by it's x/y indices, find all neighbors in towers that are given by an area
   * (min/max x/y neighbor indices).
   *
   * With the useNewton3 parameter, the lists can be either built containing all, or only the forward neighbors.
   * If an cluster A interacts with cluster B, then this interaction will either show up only once in the
   * interaction lists of the custers (for newton3 == true) or show up in the interaction lists of both
   * (for newton3 == false)
   *
   * @tparam FunType type of function
   * @param towerIndexX The x index of the given tower.
   * @param towerIndexY The y index of the given tower.
   * @param minNeighborIndexX The minimum neighbor tower index in x direction.
   * @param maxNeighborIndexX The maximum neighbor tower index in x direction.
   * @param minNeighborIndexY The minimum neighbor tower index in y direction.
   * @param maxNeighborIndexY The maximum neighbor tower index in y direction.
   * @param useNewton3 Specifies, whether neighbor lists should contain only forward neighbors.
   * @param function Function to apply on every neighbor tower. Typically this is calculateNeighborsBetweenTowers().
   */
  template <class FunType>
  void iterateNeighborTowers(const int towerIndexX, const int towerIndexY, const int minNeighborIndexX,
                             const int maxNeighborIndexX, const int minNeighborIndexY, const int maxNeighborIndexY,
                             const bool useNewton3, FunType function) {
    auto &tower = getTower(towerIndexX, towerIndexY);
    // for all neighbor towers
    for (int neighborIndexY = minNeighborIndexY; neighborIndexY <= maxNeighborIndexY; neighborIndexY++) {
      double distBetweenTowersY = std::max(0, std::abs(towerIndexY - neighborIndexY) - 1) * _towerSideLength;

      for (int neighborIndexX = minNeighborIndexX; neighborIndexX <= maxNeighborIndexX; neighborIndexX++) {
        if (useNewton3 and not isForwardNeighbor(towerIndexX, towerIndexY, neighborIndexX, neighborIndexY)) {
          continue;
        }

        double distBetweenTowersX = std::max(0, std::abs(towerIndexX - neighborIndexX) - 1) * _towerSideLength;

        // calculate distance in xy-plane
        auto distBetweenTowersXYsqr = distBetweenTowersX * distBetweenTowersX + distBetweenTowersY * distBetweenTowersY;
        // skip if already longer than interactionLength
        if (distBetweenTowersXYsqr <= _interactionLengthSqr) {
          auto &neighborTower = getTower(neighborIndexX, neighborIndexY);

          function(tower, neighborTower, distBetweenTowersXYsqr, useNewton3);
        }
      }
    }
  }

  /**
   * Returns the index of a imagined interaction cell with side length equal the interaction length that contains the
   * given tower.
   * @param towerIndexX The x index of the given tower.
   * @param towerIndexY The y index of the given tower.
   * @return The index of the interaction cell containing the given tower.
   */
  int get1DInteractionCellIndexForTower(const int towerIndexX, const int towerIndexY) {
    const int interactionCellTowerX = towerIndexX / _interactionLengthInTowers;
    const int interactionCellTowerY = towerIndexY / _interactionLengthInTowers;

    const int numInteractionCellsX = static_cast<int>(std::ceil(_towersPerDim[0] / (double)_interactionLengthInTowers));

    return interactionCellTowerX + numInteractionCellsX * interactionCellTowerY;
  }

  /**
   * Decides if a given neighbor tower is a forward neighbor to a given tower.
   * A forward neighbor is either in a interaction cell with a higher index
   * or in the same interaction cell with a higher tower index.
   *
   * Helps the VCLC06Traversal to have no data races.
   *
   * @param towerIndexX The x-index of the given tower.
   * @param towerIndexY The y-index of the given tower.
   * @param neighborIndexX The x-index of the given neighbor tower.
   * @param neighborIndexY The y-index of the given neighbor tower.
   * @return True, if neighbor is a forward neighbor of tower.
   */
  bool isForwardNeighbor(const int towerIndexX, const int towerIndexY, const int neighborIndexX,
                         const int neighborIndexY) {
    auto interactionCellTowerIndex1D = get1DInteractionCellIndexForTower(towerIndexX, towerIndexY);
    auto interactionCellNeighborIndex1D = get1DInteractionCellIndexForTower(neighborIndexX, neighborIndexY);

    if (interactionCellNeighborIndex1D > interactionCellTowerIndex1D) {
      return true;
    } else if (interactionCellNeighborIndex1D < interactionCellTowerIndex1D) {
      return false;
    }  // else if (interactionCellNeighborIndex1D == interactionCellTowerIndex1D) ...

    auto towerIndex1D = towerIndex2DTo1D(towerIndexX, towerIndexY);
    auto neighborIndex1D = towerIndex2DTo1D(neighborIndexX, neighborIndexY);

    return neighborIndex1D >= towerIndex1D;
  }

  /**
   * Calculates for all clusters in the given tower:
   *    - all neighbor clusters within the interaction length that are contained in the given neighbor tower.
   *
   * @param towerA The given tower.
   * @param towerB The given neighbor tower.
   * @param distBetweenTowersXYsqr The distance in the xy-plane between the towers.
   * @param useNewton3 Specifies, whether neighbor lists should use newton3. This changes the way what the lists
   * contain. If an cluster A interacts with cluster B, then this interaction will either show up only once in the
   * interaction lists of the custers (for newton3 == true) or show up in the interaction lists of both (for newton3 ==
   * false)
   */
  void calculateNeighborsBetweenTowers(internal::ClusterTower<Particle> &towerA,
                                       internal::ClusterTower<Particle> &towerB, double distBetweenTowersXYsqr,
                                       bool useNewton3) {
    const bool isSameTower = (&towerA == &towerB);
    for (size_t clusterIndexInTowerA = 0; clusterIndexInTowerA < towerA.getNumClusters(); clusterIndexInTowerA++) {
      // if we are within one tower depending on newton3 only look at forward neighbors
      auto startClusterIndexInTowerB = isSameTower and useNewton3 ? clusterIndexInTowerA + 1 : 0;
      auto &clusterA = towerA.getCluster(clusterIndexInTowerA);
      auto [clusterABoxBottom, clusterABoxTop, clusterAContainsParticles] = clusterA.getZMinMax();

      if (clusterAContainsParticles) {
        for (size_t clusterIndexInTowerB = startClusterIndexInTowerB; clusterIndexInTowerB < towerB.getNumClusters();
             clusterIndexInTowerB++) {
          // a cluster cannot be a neighbor to itself
          // If newton3 is true this is not possible because of the choice of the start index.
          if (not useNewton3 and isSameTower and clusterIndexInTowerA == clusterIndexInTowerB) {
            continue;
          }
          auto &clusterB = towerB.getCluster(clusterIndexInTowerB);
          auto [clusterBBoxBottom, clusterBBoxTop, clusterBcontainsParticles] = clusterB.getZMinMax();
          if (clusterBcontainsParticles) {
            double distZ = bboxDistance(clusterABoxBottom, clusterABoxTop, clusterBBoxBottom, clusterBBoxTop);
            if (distBetweenTowersXYsqr + distZ * distZ <= _interactionLengthSqr) {
              clusterA.addNeighbor(clusterB);
            }
          }
        }
      }
    }
  }

  /**
   * Calculates the distance of two bounding boxes in one dimension. Assumes disjoint bounding boxes.
   *
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

  /**
   * Returns the tower the given 3D coordinates are in.
   * If the location is outside of the domain, the tower nearest tower is returned.
   *
   * @param location The 3D coordinates.
   * @return Tower reference.
   */
  auto &getTower(std::array<double, 3> location) {
    auto [towerIndexX, towerIndexY] = getTowerCoordinates(location);
    return getTower(towerIndexX, towerIndexY);
  }

  /**
   * Returns the coordinates of the tower in the tower grid the given 3D coordinates are in.
   * If the location is outside of the domain, the tower nearest tower is returned.
   *
   * @param location The 3D coordinates.
   * @return Tower reference.
   */
  std::array<size_t, 2> getTowerCoordinates(std::array<double, 3> location) {
    std::array<size_t, 2> towerIndex2D{};

    for (int dim = 0; dim < 2; dim++) {
      const auto towerDimIndex =
          static_cast<long int>(floor((location[dim] - _haloBoxMin[dim]) * _towerSideLengthReciprocal));
      const auto towerDimIndexNonNegative = static_cast<size_t>(std::max(towerDimIndex, 0l));
      const auto towerDimIndexNonLargerValue = std::min(towerDimIndexNonNegative, _towersPerDim[dim] - 1);
      towerIndex2D[dim] = towerDimIndexNonLargerValue;
      /// @todo this is a sanity check to prevent doubling of particles, but could be done better! e.g. by border and
      // flag manager
      if (location[dim] >= _haloBoxMax[dim]) {
        towerIndex2D[dim] = _towersPerDim[dim] - 1;
      } else if (location[dim] < _haloBoxMin[dim]) {
        towerIndex2D[dim] = 0;
      }
    }

    return towerIndex2D;
  }

  /**
   * Returns the 1D index for the given 2D tower index.
   *
   * @param x The x-index of the tower.
   * @param y The y-index of the tower.
   * @return 1D index for _towers vector.
   */
  size_t towerIndex2DTo1D(const size_t x, const size_t y) {
    // It is necessary to use the static method in VerletClusterLists here instead of the member method, because
    // _towersPerDim does not have the new value yet in the container.
    return VerletClusterLists<Particle>::towerIndex2DTo1D(x, y, _towersPerDim);
  }

  /**
   * Returns the tower for the given 2D tower index.
   * @param x The x-index of the tower.
   * @param y The y-index of the tower.
   * @return Tower reference.
   */
  auto &getTower(const size_t x, const size_t y) { return _towers[towerIndex2DTo1D(x, y)]; }
};

}  //  namespace internal

}  //  namespace autopas
