/**
 * @file VerletClusterListsRebuilder.h
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
 *
 * @note Towers are always built on the xy plane towering into the z dimension.
 *
 * @tparam Particle The type of the particle the container contains.
 */
template <class Particle>
class VerletClusterListsRebuilder {
 public:
  /**
   * Type alias for the neighbor list buffer.
   */
  using NeighborListsBuffer_T = NeighborListsBuffer<const internal::Cluster<Particle> *, internal::Cluster<Particle> *>;

 private:
  size_t _clusterSize;

  NeighborListsBuffer_T &_neighborListsBuffer;

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
   * @param neighborListsBuffer Buffer structure to hold all neighbor lists.
   * @param clusterSize Size of the clusters in particles.
   */
  VerletClusterListsRebuilder(const VerletClusterLists<Particle> &clusterList,
                              std::vector<ClusterTower<Particle>> &towers, std::vector<Particle> &particlesToAdd,
                              NeighborListsBuffer_T &neighborListsBuffer, size_t clusterSize)
      : _clusterSize(clusterSize),
        _neighborListsBuffer(neighborListsBuffer),
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
    using namespace autopas::utils::ArrayMath::literals;
    // get rid of dummies
    for (auto &tower : _towers) {
      tower.deleteDummyParticles();
    }

    // count particles by accumulating tower sizes
    const size_t numParticles =
        std::accumulate(_towers.begin(), _towers.end(), _particlesToAdd.size(), [](auto acc, const auto &tower) {
          // actually we want only the number of owned or halo particles but dummies were just deleted.
          return acc + tower.getNumAllParticles();
        });

    // calculate new number of towers and their size
    const auto boxSizeWithHalo = _haloBoxMax - _haloBoxMin;
    _towerSideLength = estimateOptimalGridSideLength(numParticles, boxSizeWithHalo, _clusterSize);
    _towerSideLengthReciprocal = 1 / _towerSideLength;
    _interactionLengthInTowers = static_cast<int>(std::ceil(_interactionLength * _towerSideLengthReciprocal));
    _towersPerDim = calculateTowersPerDim(boxSizeWithHalo, _towerSideLengthReciprocal);
    const size_t numTowersNew = _towersPerDim[0] * _towersPerDim[1];

    // collect all particles that are now not in the right tower anymore
    auto invalidParticles = collectOutOfBoundsParticlesFromTowers();
    // collect all remaining particles that are not yet assigned to towers
    invalidParticles.push_back(std::move(_particlesToAdd));
    _particlesToAdd.clear();
    const auto numTowersOld = _towers.size();
    // if we have less towers than before, collect all particles from the unused towers.
    for (size_t i = numTowersNew; i < numTowersOld; ++i) {
      invalidParticles.push_back(std::move(_towers[i].particleVector()));
    }

    // resize to number of towers.
    // Attention! This uses the dummy constructor so we still need to set the desired cluster size.
    _towers.resize(numTowersNew);

    // create more towers if needed and make an estimate for how many particles memory needs to be allocated
    // Factor is more or less a random guess.
    // Historically 2.7 used to be good but in other tests there was no significant difference to lower values.
    const auto sizeEstimation =
        static_cast<size_t>((static_cast<double>(numParticles) / static_cast<double>(numTowersNew)) * 1.2);
    for (auto &tower : _towers) {
      // Set potentially new towers to the desired cluster size
      tower.setClusterSize(_clusterSize);
      tower.reserve(sizeEstimation);
    }

    sortParticlesIntoTowers(invalidParticles);

    // estimate the number of clusters by particles divided by cluster size + one extra per tower.
    _neighborListsBuffer.reserveNeighborLists(numParticles / _clusterSize + numTowersNew);
    // generate clusters and count them
    size_t numClusters = 0;
    for (auto &tower : _towers) {
      numClusters += tower.generateClusters();
      for (auto &cluster : tower.getClusters()) {
        // VCL stores the references to the lists in the clusters, therefore there is no need to create a
        // cluster -> list lookup structure in the buffer structure
        const auto listID = _neighborListsBuffer.addNeighborList();
        cluster.setNeighborList(&(_neighborListsBuffer.template getNeighborListRef<false>(listID)));
      }
    }

    return std::make_tuple(_towerSideLength, _interactionLengthInTowers, _towersPerDim, numClusters);
  }

  /**
   * Rebuilds the neighbor lists and fills Clusters with dummies as described in
   * ClusterTower::setDummyValues.
   * @param useNewton3 Specifies, whether neighbor lists should use newton3. This changes the way what the lists
   * contain. If an cluster A interacts with cluster B, then this interaction will either show up only once in the
   * interaction lists of the custers (for newton3 == true) or show up in the interaction lists of both (for newton3 ==
   * false)
   */
  void rebuildNeighborListsAndFillClusters(bool useNewton3) {
    clearNeighborListsAndmoveDummiesIntoClusters();
    updateNeighborLists(useNewton3);

    double dummyParticleDistance = _interactionLength * 2;
    double startDummiesX = 1000 * _haloBoxMax[0];
    for (size_t index = 0; index < _towers.size(); index++) {
      _towers[index].setDummyValues(startDummiesX + static_cast<double>(index) * dummyParticleDistance,
                                    dummyParticleDistance);
    }
  }
  /**
   * Estimates the optimal grid side length.
   * @param numParticles The number of particles in the container.
   * @param boxSize The size of the domain.
   * @param clusterSize the number of particles per cluster.
   * @return an estimated optimal grid side length.
   */
  [[nodiscard]] static double estimateOptimalGridSideLength(size_t numParticles, const std::array<double, 3> &boxSize,
                                                            size_t clusterSize) {
    const double volume = boxSize[0] * boxSize[1] * boxSize[2];
    if (numParticles > 0) {
      // estimate particle density
      const double density = static_cast<double>(numParticles) / volume;

      return std::cbrt(static_cast<double>(clusterSize) / density);
    } else {
      return std::max(boxSize[0], boxSize[1]);
    }
  }

  /**
   * Calculates the cells per dimension in the container using the _towerSideLengthReciprocal.
   * @param boxSize the size of the domain.
   * @param towerSideLengthReciprocal 1.0 / towerSidelength.
   * @return the cells per dimension in the container.
   */
  [[nodiscard]] static std::array<size_t, 2> calculateTowersPerDim(const std::array<double, 3> &boxSize,
                                                                   double towerSideLengthReciprocal) {
    std::array<size_t, 2> towersPerDim{};
    for (int d = 0; d < 2; d++) {
      towersPerDim[d] = static_cast<size_t>(std::ceil(boxSize[d] * towerSideLengthReciprocal));
      // at least one cell
      towersPerDim[d] = std::max(towersPerDim[d], 1ul);
    }
    return towersPerDim;
  }

  /**
   * Clears previously saved neighbors from clusters and sets the 3D positions of the dummy particles to inside of the
   * cluster to avoid all dummies being in one place and potentially trigger cluster-cluster distance evaluations.
   */
  void clearNeighborListsAndmoveDummiesIntoClusters() {
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
   * Collect all particles that are stored in the wrong towers.
   * The particles are deleted from their towers.
   *
   * @return
   */
  std::vector<std::vector<Particle>> collectOutOfBoundsParticlesFromTowers() {
    std::vector<std::vector<Particle>> outOfBoundsParticles;
    outOfBoundsParticles.resize(_towers.size());
    for (size_t towerIndex = 0; towerIndex < _towers.size(); towerIndex++) {
      const auto &[towerBoxMin, towerBoxMax] = VerletClusterLists<Particle>::getTowerBoundingBox(towerIndex);
      outOfBoundsParticles[towerIndex] = _towers[towerIndex].collectOutOfBoundsParticles(towerBoxMin, towerBoxMax);
    }
    return outOfBoundsParticles;
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
    /// @todo: find sensible chunk size
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t index = 0; index < numVectors; index++) {
      const std::vector<Particle> &vector = particles2D[index];
      for (const auto &particle : vector) {
        if (utils::inBox(particle.getR(), _haloBoxMin, _haloBoxMax)) {
          auto &tower = getTower(particle.getR());
          tower.addParticle(particle);
        } else {
          AutoPasLog(TRACE, "Not adding particle to VerletClusterLists container, because it is outside the halo:\n{}",
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

    const int numInteractionCellsX =
        static_cast<int>(std::ceil(_towersPerDim[0] / static_cast<double>(_interactionLengthInTowers)));

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
    const auto interactionCellTowerIndex1D = get1DInteractionCellIndexForTower(towerIndexX, towerIndexY);
    const auto interactionCellNeighborIndex1D = get1DInteractionCellIndexForTower(neighborIndexX, neighborIndexY);

    if (interactionCellNeighborIndex1D > interactionCellTowerIndex1D) {
      return true;
    } else if (interactionCellNeighborIndex1D < interactionCellTowerIndex1D) {
      return false;
    }  // else if (interactionCellNeighborIndex1D == interactionCellTowerIndex1D) ...

    const auto towerIndex1D = towerIndex2DTo1D(towerIndexX, towerIndexY);
    const auto neighborIndex1D = towerIndex2DTo1D(neighborIndexX, neighborIndexY);

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
    const auto interactionLengthFracOfDomainZ = _interactionLength / (_haloBoxMax[2] - _haloBoxMin[2]);
    // Seems to find a good middle ground between not too much memory allocated and no additional allocations
    // when calling clusterA.addNeighbor(clusterB)
    const auto neighborListReserveHeuristicFactor = (interactionLengthFracOfDomainZ * 2.1) / _clusterSize;
    const bool isSameTower = (&towerA == &towerB);
    for (size_t clusterIndexInTowerA = 0; clusterIndexInTowerA < towerA.getNumClusters(); clusterIndexInTowerA++) {
      // if we are within one tower depending on newton3 only look at forward neighbors
      const auto startClusterIndexInTowerB = isSameTower and useNewton3 ? clusterIndexInTowerA + 1ul : 0ul;
      auto &clusterA = towerA.getCluster(clusterIndexInTowerA);
      const auto [clusterABoxBottom, clusterABoxTop, clusterAContainsParticles] = clusterA.getZMinMax();

      if (clusterAContainsParticles) {
        clusterA.getNeighbors()->reserve((towerA.numParticles() + 8 * towerB.numParticles()) *
                                         neighborListReserveHeuristicFactor);
        for (size_t clusterIndexInTowerB = startClusterIndexInTowerB; clusterIndexInTowerB < towerB.getNumClusters();
             clusterIndexInTowerB++) {
          // a cluster cannot be a neighbor to itself
          // If newton3 is true this is not possible because of the choice of the start index.
          if (not useNewton3 and isSameTower and clusterIndexInTowerA == clusterIndexInTowerB) {
            continue;
          }
          // can't be const because it will potentially be added as a non-const neighbor
          auto &clusterB = towerB.getCluster(clusterIndexInTowerB);
          const auto [clusterBBoxBottom, clusterBBoxTop, clusterBcontainsParticles] = clusterB.getZMinMax();
          if (clusterBcontainsParticles) {
            const double distZ = bboxDistance(clusterABoxBottom, clusterABoxTop, clusterBBoxBottom, clusterBBoxTop);
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
  auto &getTower(const std::array<double, 3> &location) {
    auto [towerIndexX, towerIndexY] = getTowerCoordinates(location);
    return getTower(towerIndexX, towerIndexY);
  }

  /**
   * Returns the 2D index of the tower in the tower grid the given 3D coordinates are in.
   * If the location is outside of the domain, the tower nearest tower is returned.
   *
   * @param location The 3D coordinates.
   * @return 2D tower index.
   */
  [[nodiscard]] std::array<size_t, 2> getTowerCoordinates(const std::array<double, 3> &location) const {
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
  [[nodiscard]] size_t towerIndex2DTo1D(const size_t x, const size_t y) const {
    // It is necessary to use the static method in VerletClusterLists here instead of the member method, because
    // _towersPerDim does not have the new value yet in the container.
    return VerletClusterLists<Particle>::towerIndex2DTo1D(x, y, _towersPerDim);
  }

  /**
   * Returns the 2D index for the given 1D index of a tower. Static version.
   *
   * @param index
   * @return the 2D index for the given 1D index of a tower.
   */
  [[nodiscard]] std::array<size_t, 2> towerIndex1DTo2D(const size_t index) const {
    // It is necessary to use the static method in VerletClusterLists here instead of the member method, because
    // _towersPerDim does not have the new value yet in the container.
    return VerletClusterLists<Particle>::towerIndex1DTo2D(index, _towersPerDim[0]);
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
