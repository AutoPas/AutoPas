/**
 * @file VerletClusterListsRebuilder.h
 * @author humig
 * @date 29.07.19
 */

#pragma once

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/inBox.h"

namespace autopas::internal {

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
  ClusterTowerBlock2D<Particle> &_towerBlock;
  double _interactionLengthSqr;
  bool _newton3;

 public:
  /**
   * Constructs the builder from the cluster list.
   *
   * @param towerBlock The towers from the cluster list to rebuild.
   * @param particlesToAdd New particles to add.
   * @param neighborListsBuffer Buffer structure to hold all neighbor lists.
   * @param clusterSize Size of the clusters in particles.
   * @param interactionLengthSqr Squared interaction length (cutoff + skin)^2
   * @param newton3 If the current configuration uses newton3
   */
  VerletClusterListsRebuilder(ClusterTowerBlock2D<Particle> &towerBlock, std::vector<Particle> &particlesToAdd,
                              NeighborListsBuffer_T &neighborListsBuffer, size_t clusterSize,
                              double interactionLengthSqr, bool newton3)
      : _clusterSize(clusterSize),
        _neighborListsBuffer(neighborListsBuffer),
        _particlesToAdd(particlesToAdd),
        _towerBlock(towerBlock),
        _interactionLengthSqr(interactionLengthSqr),
        _newton3(newton3) {}

  /**
   * Rebuilds the towers, sorts the particles into them and creates the clusters with a reference to an uninitialized
   * neighbor list.
   *
   * @return number of clusters
   */
  size_t rebuildTowersAndClusters() {
    using namespace autopas::utils::ArrayMath::literals;
    // get rid of dummies
    for (auto &tower : _towerBlock) {
      tower.deleteDummyParticles();
    }

    // count particles by accumulating tower sizes
    const size_t numParticles = std::accumulate(_towerBlock.begin(), _towerBlock.end(), _particlesToAdd.size(),
                                                [](auto acc, const auto &tower) {
                                                  // actually we want only the number of owned or halo particles
                                                  // but dummies were just deleted.
                                                  return acc + tower.getNumActualParticles();
                                                });

    // calculate new number of towers and their size
    const auto boxSizeWithHalo = _towerBlock.getHaloBoxMax() - _towerBlock.getHaloBoxMin();
    const auto numTowersOld = _towerBlock.size();
    const auto [towerSideLength, numTowersPerDim] =
        _towerBlock.estimateOptimalGridSideLength(numParticles, _clusterSize);
    const auto numTowersNew = numTowersPerDim[0] * numTowersPerDim[1];

    // collect all particles that are now not in the right tower anymore
    auto invalidParticles = collectOutOfBoundsParticlesFromTowers();
    // collect all remaining particles that are not yet assigned to towers
    invalidParticles.push_back(std::move(_particlesToAdd));
    _particlesToAdd.clear();
    // if we have less towers than before, collect all particles from the unused towers.
    for (size_t i = numTowersNew; i < numTowersOld; ++i) {
      invalidParticles.push_back(std::move(_towerBlock[i].particleVector()));
    }

    // resize to number of towers.
    // Attention! This uses the dummy constructor so we still need to set the desired cluster size.
    _towerBlock.resize(towerSideLength, numTowersPerDim);

    // create more towers if needed and make an estimate for how many particles memory needs to be allocated
    // Factor is more or less a random guess.
    // Historically 2.7 used to be good but in other tests there was no significant difference to lower values.
    const auto sizeEstimation =
        static_cast<size_t>((static_cast<double>(numParticles) / static_cast<double>(numTowersNew)) * 1.2);
    for (auto &tower : _towerBlock) {
      // Set potentially new towers to the desired cluster size
      tower.setClusterSize(_clusterSize);
      tower.reserve(sizeEstimation);
    }

    sortParticlesIntoTowers(invalidParticles);

    // estimate the number of clusters by particles divided by cluster size + one extra per tower.
    _neighborListsBuffer.reserveNeighborLists(numParticles / _clusterSize + numTowersNew);
    // generate clusters and count them
    size_t numClusters = 0;
    for (auto &tower : _towerBlock) {
      numClusters += tower.generateClusters();
      for (auto clusterIter = _newton3 ? tower.getClusters().begin() : tower.getFirstOwnedCluster();
           clusterIter < (_newton3 ? tower.getClusters().end() : tower.getFirstTailHaloCluster()); ++clusterIter) {
        // VCL stores the references to the lists in the clusters, therefore there is no need to create a
        // cluster -> list lookup structure in the buffer structure
        const auto listID = _neighborListsBuffer.getNewNeighborList();
        clusterIter->setNeighborList(&(_neighborListsBuffer.template getNeighborListRef<false>(listID)));
      }
    }

    return numClusters;
  }

  /**
   * Rebuilds the neighbor lists and fills Clusters with dummies as described in ClusterTower::setDummyValues.
   *
   * @note Here, _newton3 decides, whether neighbor lists should use newton3. This changes what the lists contain.
   * For two interacting clusters A and B, if _newton3 == false, the interaction A->B is in the list of cluster B, and
   * B->A is in cluster A. If _newton3 == true, the two-way interaction A<->B will only be in the cluster that comes
   * first when iterating through towers.
   */
  void rebuildNeighborListsAndFillClusters() {
    clearNeighborListsAndMoveDummiesIntoClusters();
    updateNeighborLists();

    double dummyParticleDistance = _towerBlock.getInteractionLength() * 2;
    double startDummiesX = 1000 * _towerBlock.getHaloBoxMax()[0];
    for (size_t index = 0; index < _towerBlock.size(); index++) {
      _towerBlock[index].setDummyValues(startDummiesX + static_cast<double>(index) * dummyParticleDistance,
                                        dummyParticleDistance);
    }
  }

  /**
   * Clears previously saved neighbors from clusters and sets the 3D positions of the dummy particles to inside of the
   * cluster to avoid all dummies being in one place and potentially trigger cluster-cluster distance evaluations.
   */
  void clearNeighborListsAndMoveDummiesIntoClusters() {
    for (auto &tower : _towerBlock) {
      tower.setDummyParticlesToLastActualParticle();
      for (auto clusterIter = tower.getFirstOwnedCluster(); clusterIter < tower.getFirstTailHaloCluster();
           ++clusterIter) {
        clusterIter->clearNeighbors();
      }
    }
  }

  /**
   * Takes all particles from all towers and returns them. Towers are cleared afterwards.
   * @return All particles in the container sorted in 2D as they were in the towers.
   */
  std::vector<std::vector<Particle>> collectAllParticlesFromTowers() {
    std::vector<std::vector<Particle>> invalidParticles;
    invalidParticles.resize(_towerBlock.size());
    for (size_t towerIndex = 0; towerIndex < _towerBlock.size(); towerIndex++) {
      invalidParticles[towerIndex] = _towerBlock[towerIndex].collectAllActualParticles();
      _towerBlock[towerIndex].clear();
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
    outOfBoundsParticles.resize(_towerBlock.size());
    for (size_t towerIndex = 0; towerIndex < _towerBlock.size(); towerIndex++) {
      const auto &[towerBoxMin, towerBoxMax] = _towerBlock.getTowerBoundingBox(towerIndex);
      outOfBoundsParticles[towerIndex] = _towerBlock[towerIndex].collectOutOfBoundsParticles(towerBoxMin, towerBoxMax);
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
        if (utils::inBox(particle.getR(), _towerBlock.getHaloBoxMin(), _towerBlock.getHaloBoxMax())) {
          auto &tower = _towerBlock.getTowerAtPosition(particle.getR());
          tower.addParticle(particle);
        } else {
          AutoPasLog(TRACE, "Not adding particle to VerletClusterLists container, because it is outside the halo:\n{}",
                     particle.toString());
        }
      }
    }
  }

  /**
   * Updates the neighbor lists for all clusters.
   */
  void updateNeighborLists() {
    const int maxTowerIndexX = _towerBlock.getTowersPerDim()[0] - 1;
    const int maxTowerIndexY = _towerBlock.getTowersPerDim()[1] - 1;
    const auto numTowersPerInteractionLength = _towerBlock.getNumTowersPerInteractionLength();
    // for all towers
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (int towerIndexY = 0; towerIndexY <= maxTowerIndexY; towerIndexY++) {
      for (int towerIndexX = 0; towerIndexX <= maxTowerIndexX; towerIndexX++) {
        // calculate extent of interesting tower 2D indices
        const int minX = std::max(towerIndexX - numTowersPerInteractionLength, 0);
        const int minY = std::max(towerIndexY - numTowersPerInteractionLength, 0);
        const int maxX = std::min(towerIndexX + numTowersPerInteractionLength, maxTowerIndexX);
        const int maxY = std::min(towerIndexY + numTowersPerInteractionLength, maxTowerIndexY);

        iterateNeighborTowers(towerIndexX, towerIndexY, minX, maxX, minY, maxY,
                              [this](auto &towerA, auto &towerB) { calculateNeighborsBetweenTowers(towerA, towerB); });
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
   * @note _newton3 Specifies, whether neighbor lists should contain only forward neighbors.
   *
   * @tparam FunType type of function
   * @param towerIndexX The x index of the given tower.
   * @param towerIndexY The y index of the given tower.
   * @param minNeighborIndexX The minimum neighbor tower index in x direction.
   * @param maxNeighborIndexX The maximum neighbor tower index in x direction.
   * @param minNeighborIndexY The minimum neighbor tower index in y direction.
   * @param maxNeighborIndexY The maximum neighbor tower index in y direction.
   * @param function Function to apply on every neighbor tower. Typically this is calculateNeighborsBetweenTowers().
   */
  template <class FunType>
  void iterateNeighborTowers(const int towerIndexX, const int towerIndexY, const int minNeighborIndexX,
                             const int maxNeighborIndexX, const int minNeighborIndexY, const int maxNeighborIndexY,
                             FunType function) {
    auto &tower = _towerBlock.getTowerByIndex2D(towerIndexX, towerIndexY);
    // for all neighbor towers
    for (int neighborIndexY = minNeighborIndexY; neighborIndexY <= maxNeighborIndexY; neighborIndexY++) {
      double distBetweenTowersY =
          std::max(0, std::abs(towerIndexY - neighborIndexY) - 1) * _towerBlock.getTowerSideLength()[1];

      for (int neighborIndexX = minNeighborIndexX; neighborIndexX <= maxNeighborIndexX; neighborIndexX++) {
        if (_newton3 and not isForwardNeighbor(towerIndexX, towerIndexY, neighborIndexX, neighborIndexY)) {
          continue;
        }

        double distBetweenTowersX =
            std::max(0, std::abs(towerIndexX - neighborIndexX) - 1) * _towerBlock.getTowerSideLength()[0];

        // calculate distance in xy-plane
        auto distBetweenTowersXYsqr = distBetweenTowersX * distBetweenTowersX + distBetweenTowersY * distBetweenTowersY;
        // skip if already longer than interactionLength
        if (distBetweenTowersXYsqr <= _interactionLengthSqr) {
          auto &neighborTower = _towerBlock.getTowerByIndex2D(neighborIndexX, neighborIndexY);

          function(tower, neighborTower);
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
    const auto numTowersPerInteractionLength = _towerBlock.getNumTowersPerInteractionLength();
    const int interactionCellTowerX = towerIndexX / numTowersPerInteractionLength;
    const int interactionCellTowerY = towerIndexY / numTowersPerInteractionLength;

    const int numInteractionCellsX = static_cast<int>(
        std::ceil(_towerBlock.getTowersPerDim()[0] / static_cast<double>(numTowersPerInteractionLength)));

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

    const auto towerIndex1D = _towerBlock.towerIndex2DTo1D(towerIndexX, towerIndexY);
    const auto neighborIndex1D = _towerBlock.towerIndex2DTo1D(neighborIndexX, neighborIndexY);

    return neighborIndex1D >= towerIndex1D;
  }

  /**
   * Calculates for all clusters in the given tower:
   *    - all neighbor clusters within the interaction length that are contained in the given neighbor tower.
   *
   * @param towerA The given tower.
   * @param towerB The given neighbor tower.
   * contain. If an cluster A interacts with cluster B, then this interaction will either show up only once in the
   * interaction lists of the custers (for newton3 == true) or show up in the interaction lists of both (for newton3 ==
   * false)
   */
  void calculateNeighborsBetweenTowers(internal::ClusterTower<Particle> &towerA,
                                       internal::ClusterTower<Particle> &towerB) {
    using autopas::utils::ArrayMath::boxDistanceSquared;

    const auto interactionLengthFracOfDomainZ =
        _towerBlock.getInteractionLength() / (_towerBlock.getHaloBoxMax()[2] - _towerBlock.getHaloBoxMin()[2]);
    const bool isSameTower = (&towerA == &towerB);
    // This heuristic seems to find a good middle ground between not too much memory allocated and no additional
    // allocations when calling clusterA.addNeighbor(clusterB)
    const auto neighborListReserveHeuristicFactor = (interactionLengthFracOfDomainZ * 2.1) / _clusterSize;
    const auto isHaloCluster = [](const auto &clusterIter, const auto &tower) {
      return clusterIter < tower.getFirstOwnedCluster() or clusterIter >= tower.getFirstTailHaloCluster();
    };
    // iterate over all clusters from tower A. In newton3 mode go over all of them, otherwise only owned.
    for (auto clusterIterA = _newton3 ? towerA.getClusters().begin() : towerA.getFirstOwnedCluster();
         clusterIterA < (_newton3 ? towerA.getClusters().end() : towerA.getFirstTailHaloCluster()); ++clusterIterA) {
      if (not clusterIterA->empty()) {
        clusterIterA->getNeighbors()->reserve((towerA.getNumActualParticles() + 8 * towerB.getNumActualParticles()) *
                                              neighborListReserveHeuristicFactor);

        // if we are within one tower depending on newton3 only look at forward neighbors
        // clusterIterB can't be const because it will potentially be added as a non-const neighbor
        for (auto clusterIterB = isSameTower and _newton3 ? clusterIterA + 1 : towerB.getClusters().begin();
             clusterIterB < towerB.getClusters().end(); ++clusterIterB) {
          // a cluster cannot be a neighbor to itself
          if (&*clusterIterA == &*clusterIterB) {
            continue;
          }
          // never do halo-halo interactions
          if (isHaloCluster(clusterIterA, towerA) and isHaloCluster(clusterIterB, towerB)) {
            continue;
          }
          if (not clusterIterB->empty()) {
            const auto [aMin, aMax] = clusterIterA->getBoundingBox();
            const auto [bMin, bMax] = clusterIterB->getBoundingBox();
            const auto boxDistSquared = boxDistanceSquared(aMin, aMax, bMin, bMax);
            if (boxDistSquared <= _interactionLengthSqr) {
              clusterIterA->addNeighbor(*clusterIterB);
            }
          }
        }
      }
    }
  }
  /**
   * Getter
   * @return
   */
  bool getNewton3() const { return _newton3; }
};
}  // namespace autopas::internal
