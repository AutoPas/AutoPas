/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <cmath>
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "autopas/containers/verletClusterLists/VerletClusterMaths.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/inBox.h"

namespace autopas {

template <class Particle>
class VerletClustersTraversalInterface;

/**
 * Particles are divided into clusters.
 * The VerletClusterLists class uses neighborhood lists for each cluster
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * @tparam Particle
 */
template <class Particle>
class VerletClusterLists : public ParticleContainer<Particle, FullParticleCell<Particle>> {
 public:
  /**
   * The number of particles in a full cluster. If necessary, constexpr can be removed it can be made a member variable.
   */
  static constexpr size_t clusterSize = 4;

  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly the
   * same side length.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param clusterSize size of clusters
   */
  VerletClusterLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                     double skin = 0)
      : ParticleContainer<Particle, FullParticleCell<Particle>>(boxMin, boxMax, cutoff, skin),
        _numClusters(0),
        _interactionLengthSqr(this->getInteractionLength() * this->getInteractionLength()),
        _interactionLengthInTowers(0),
        _neighborListIsNewton3(false) {}

  ContainerOption getContainerType() override { return ContainerOption::verletClusterLists; }

  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", utils::StringUtils::to_string(traversal->getTraversalType()));

    auto *traversalInterface = dynamic_cast<VerletClustersTraversalInterface<Particle> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setClusterLists(*this);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletClusterLists::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }

    utils::Timer timer;
    static double initTime = 0, traversalTime = 0, endTime = 0;
    timer.start();
    traversal->initTraversal();
    std::cout << (initTime += timer.stop()) << " for init!" << std::endl;
    timer.start();
    traversal->traverseParticlePairs();
    std::cout << (traversalTime += timer.stop()) << " for traversal!" << std::endl;
    timer.start();
    traversal->endTraversal();
    std::cout << (endTime += timer.stop()) << " for end!" << std::endl;
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle &p) override {
    // add particle somewhere, because lists will be rebuild anyways
    _towers[0].addParticle(p);
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(Particle &haloParticle) override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.addHaloParticle not yet implemented.");
  }

  bool updateHaloParticle(Particle &haloParticle) override { throw std::runtime_error("not yet implemented"); }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    // quick and dirty: iterate over all particles and delete halo particles
    // @todo: make this proper
    for (auto iter = this->begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
      if (not iter->isOwned()) {
        iter.deleteCurrentParticle();
      }
    }
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<Particle> updateContainer() override {
    // TODO What happens when some particles are just deleted here?
    AutoPasLog(debug, "updating container");
    // first delete all particles
    this->deleteHaloParticles();

    // next find invalid particles
    std::vector<Particle> invalidParticles;
    // @todo: parallelize
    for (auto iter = this->begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      if (not utils::inBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        invalidParticles.push_back(*iter);
        iter.deleteCurrentParticle();
      }
    }

    return invalidParticles;
  }

  bool isContainerUpdateNeeded() override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.isContainerUpdateNeeded not yet implemented");
    return false;
  }

  TraversalSelectorInfo getTraversalSelectorInfo() override { return TraversalSelectorInfo(_towersPerDim); }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, internal::ClusterTower<Particle, 4>>(&(this->_towers)));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle>();
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    utils::Timer timer;
    static double rebuildTime = 0;
    timer.start();
    rebuild(traversal->getUseNewton3());
    std::cout << (rebuildTime += timer.stop()) << " for rebuild!" << std::endl;
  }

  /**
   * Helper method to iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @tparam inParallel If the iteration should be executed in parallel or sequential.  See traverseClustersParallel()
   * for thread safety.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <bool inParallel, class LoopBody>
  void traverseClusters(LoopBody &&loopBody) {
    if (inParallel) {
      traverseClustersParallel<LoopBody>(std::forward<LoopBody>(loopBody));
    } else {
      traverseClustersSequential<LoopBody>(std::forward<LoopBody>(loopBody));
    }
  }

  unsigned long getNumParticles() override {
    unsigned long sum = 0;
    for (size_t index = 0; index < _towers.size(); index++) {
      sum += _towers[index].getNumParticles();
    }
    return sum;
  }

  /**
   * Returns the number of clusters in this container.
   * @return The number of clusters in this container.
   */
  auto getNumClusters() const { return _numClusters; }

  /**
   * Returns the grid side length of the grids in the container.
   * @return the grid side length of the grids in the container.
   */
  auto getTowerSideLength() const { return _towerSideLength; }

  /**
   * Returns the number of grids per dimension on the container.
   * @return the number of grids per dimension on the container.
   */
  auto getTowersPerDimension() const { return _towersPerDim; }

  /**
   * Returns the 2D grid for the XY-plane of this container that defines the cluster towers.
   * @return the grids of this container for usage in traversals.
   */
  auto &getTowers() { return _towers; }

  /**
   * Returns the number of particles in each cluster.
   * @return the number of particles in each cluster.
   */
  constexpr auto getClusterSize() const { return clusterSize; }

  /**
   * Returns the interaction length in towers. That is how many towers fit into one interaction length rounded up.
   * @return the interaciton length in towers.
   */
  auto getInteractionLengthInTowers() const { return _interactionLengthInTowers; }

 protected:
  /**
   * Helper method to sequentially iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <class LoopBody>
  void traverseClustersSequential(LoopBody &&loopBody) {
    for (size_t x = 0; x < _towersPerDim[0]; x++) {
      for (size_t y = 0; y < _towersPerDim[1]; y++) {
        size_t index = VerletClusterMaths::index1D(x, y, _towersPerDim);
        auto &tower = _towers[index];

        for (auto &cluster : tower.getClusters()) {
          loopBody(cluster);
        }
      }
    }
  }

  /**
   * Helper method to iterate over all clusters in parallel.
   *
   * It is always safe to modify the particles in the cluster that is passed to the given loop body. However, when
   * modifying particles from other clusters, the caller has to make sure that no data races occur. Particles must not
   * be added or removed during the traversal.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <class LoopBody>
  void traverseClustersParallel(LoopBody &&loopBody) {
#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (size_t x = 0; x < _towersPerDim[0]; x++) {
      for (size_t y = 0; y < _towersPerDim[1]; y++) {
        size_t index = VerletClusterMaths::index1D(x, y, _towersPerDim);
        auto &tower = _towers[index];

        for (auto &cluster : tower.getClusters()) {
          loopBody(cluster);
        }
      }
    }
  }

  /**
   * Recalculate grids and clusters, build verlet lists and pad clusters.
   * @param useNewton3 If the everything should be build using newton 3 or not.
   */
  void rebuild(bool useNewton3) {
    std::vector<Particle> invalidParticles = collectAllParticlesFromTowers();

    auto boxSize = ArrayMath::sub(this->getBoxMax(), this->getBoxMin());

    _towerSideLength = estimateOptimalGridSideLength(invalidParticles.size(), boxSize);
    _interactionLengthInTowers = static_cast<int>(std::ceil((this->getInteractionLength()) / _towerSideLength));
    _towerSideLengthReciprocal = 1 / _towerSideLength;

    _towersPerDim = calculateTowersPerDim(boxSize);
    // _towersPerDim[2] is always 1
    size_t numCells = _towersPerDim[0] * _towersPerDim[1];

    // resize to number of towers
    _towers.resize(numCells);

    sortParticlesIntoTowers(invalidParticles);

    size_t numClustersAcc = 0;
    for (auto &tower : _towers) {
      numClustersAcc += tower.generateClusters();
    }
    _numClusters = numClustersAcc;

    updateNeighborLists(useNewton3);

    // Make sure that each cluster contains a multiple of clusterSize particles, and that dummies don't interact.
    double dummyParticleDistance = this->getInteractionLength() * 2;
    double startDummiesX = 1000 * this->getBoxMax()[0];
    for (size_t index = 0; index < _towers.size(); index++) {
      _towers[index].fillUpWithDummyParticles(startDummiesX + index * dummyParticleDistance, dummyParticleDistance);
    }
  }

  /**
   * Takes all particles from all towers and returns them. Towers are cleared afterwards.
   * @return All particles in the container.
   */
  std::vector<Particle> collectAllParticlesFromTowers() {
    // TODO: Optimize so not all particles are copied, but the vectors of the tower moved here
    std::vector<Particle> invalidParticles;
    for (auto &tower : _towers) {
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
   * Calculates the cells per dimension in the container using the _gridSideLengthReciprocal.
   * @param boxSize the size of the domain.
   * @return the cells per dimension in the container.
   */
  [[nodiscard]] std::array<size_t, 3> calculateTowersPerDim(std::array<double, 3> boxSize) const {
    std::array<size_t, 3> cellsPerDim{};
    for (int d = 0; d < 2; d++) {
      cellsPerDim[d] = static_cast<size_t>(std::ceil(boxSize[d] * _towerSideLengthReciprocal));
      // at least one cell
      cellsPerDim[d] = std::max(cellsPerDim[d], 1ul);
    }
    cellsPerDim[2] = 1ul;
    return cellsPerDim;
  }

  /**
   * Sorts all passed particles in the appropriate clusters.
   * @param particles The particles to sort in the clusters.
   */
  void sortParticlesIntoTowers(std::vector<Particle> &particles) {
    for (auto &particle : particles) {
      if (utils::inBox(particle.getR(), this->getBoxMin(), this->getBoxMax())) {
        auto index = get1DIndexOfPosition(particle.getR());
        _towers[index].addParticle(particle);
      }
    }
  }

  /**
   * Updates the neighbor lists.
   *
   * @param useNewton3 If newton 3 should be used to build the neighbor lists or not. If true, only saves neighbor
   * clusters that have a higher index that the current cluster. (@see buildClusterIndexMap())
   */
  void updateNeighborLists(bool useNewton3) {
    _neighborListIsNewton3 = useNewton3;

    const int maxTowerIndexX = _towersPerDim[0] - 1;
    const int maxTowerIndexY = _towersPerDim[1] - 1;
    // for all towers
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
    auto &tower = _towers[VerletClusterMaths::index1D(towerIndexX, towerIndexY, _towersPerDim)];
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
          auto &neighborTower = _towers[VerletClusterMaths::index1D(neighborIndexX, neighborIndexY, _towersPerDim)];

          calculateNeighborsForTowerPair(tower, neighborTower, distBetweenTowersXYsqr);
        }
      }
    }
  }

  int get1DInteractionCellIndexForTower(const int towerIndexX, const int towerIndexY) {
    const int interactionCellTowerX = towerIndexX / _interactionLengthInTowers;
    const int interactionCellTowerY = towerIndexY / _interactionLengthInTowers;

    // TODO: Check if having this as member improves performance
    const int numInteractionCellsX = static_cast<int>(std::ceil(_towersPerDim[0] / (double)_interactionLengthInTowers));

    return interactionCellTowerX + numInteractionCellsX * interactionCellTowerY;
  }

  bool shouldTowerContainOtherAsNeighborWithNewton3(const int towerIndexX, const int towerIndexY,
                                                    const int neighborIndexX, const int neighborIndexY) {
    auto interactionCellTowerIndex1D = get1DInteractionCellIndexForTower(towerIndexX, towerIndexY);
    auto interactionCellNeighborIndex1D = get1DInteractionCellIndexForTower(neighborIndexX, neighborIndexY);

    return interactionCellNeighborIndex1D > interactionCellTowerIndex1D or
           (interactionCellNeighborIndex1D == interactionCellTowerIndex1D and
            VerletClusterMaths::index1D(neighborIndexX, neighborIndexY, _towersPerDim) >=
                VerletClusterMaths::index1D(towerIndexX, towerIndexY, _towersPerDim));
  }

  void calculateNeighborsForTowerPair(internal::ClusterTower<Particle, clusterSize> &tower,
                                      internal::ClusterTower<Particle, clusterSize> &neighborTower,
                                      double distBetweenTowersXYsqr) {
    for (size_t towerIndex = 0; towerIndex < tower.getNumClusters(); towerIndex++) {
      auto startIndexNeighbor = _neighborListIsNewton3 and &tower == &neighborTower ? towerIndex + 1 : 0;
      auto &towerCluster = tower.getCluster(towerIndex);
      double towerClusterBoxBottom = towerCluster.getParticle(0).getR()[2];
      double towerClusterBoxTop = towerCluster.getParticle(clusterSize - 1).getR()[2];

      for (size_t neighborIndex = startIndexNeighbor; neighborIndex < neighborTower.getNumClusters(); neighborIndex++) {
        if (not _neighborListIsNewton3 and towerIndex == neighborIndex) {
          continue;
        }
        auto &neighborCluster = neighborTower.getCluster(neighborIndex);
        double neighborClusterBoxBottom = neighborCluster.getParticle(0).getR()[2];
        double neighborClusterBoxTop = neighborCluster.getParticle(clusterSize - 1).getR()[2];

        double distZ =
            bboxDistance(towerClusterBoxBottom, towerClusterBoxTop, neighborClusterBoxBottom, neighborClusterBoxTop);
        if (distBetweenTowersXYsqr + distZ * distZ <= _interactionLengthSqr) {
          towerCluster.getNeighbors().push_back(&neighborCluster);
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

      /**
       * Gets the 1d tower index containing a particle in given position.
       * @param pos the position of the particle
       * @return the index of the tower containing the given position
       */
      [[nodiscard]] size_t get1DIndexOfPosition(const std::array<double, 3> &pos) const {
    std::array<size_t, 2> cellIndex{};

    for (int dim = 0; dim < 2; dim++) {
      const long int value =
          (static_cast<long int>(floor((pos[dim] - this->getBoxMin()[dim]) * _towerSideLengthReciprocal))) + 1l;
      const size_t nonNegativeValue = static_cast<size_t>(std::max(value, 0l));
      const size_t nonLargerValue = std::min(nonNegativeValue, _towersPerDim[dim] - 1);
      cellIndex[dim] = nonLargerValue;
      // @todo this is a sanity check to prevent doubling of particles, but could be done better! e.g. by border and
      // flag manager
      if (pos[dim] >= this->getBoxMax()[dim]) {
        cellIndex[dim] = _towersPerDim[dim] - 1;
      } else if (pos[dim] < this->getBoxMin()[dim]) {
        cellIndex[dim] = 0;
      }
    }

    return VerletClusterMaths::index1D(cellIndex[0], cellIndex[1], _towersPerDim);
  }

 private:
  /**
   * internal storage, particles are split into a grid in xy-dimension
   */
  std::vector<internal::ClusterTower<Particle, clusterSize>> _towers{1};

  /**
   * Dimensions of the 2D xy-grid.
   */
  std::array<size_t, 3> _towersPerDim{};

  /**
   * Side length of xy-grid.
   */
  double _towerSideLength{0.};

  /**
   *  Reciprocal of _gridSideLength.
   */
  double _towerSideLengthReciprocal{0.};

  /**
   * The number of clusters in the container.
   */
  size_t _numClusters;

  /**
   * (_cutoff + _skin)^2.
   */
  double _interactionLengthSqr;

  /**
   * The interaction length in number of towers it reaches.
   * static_cast<int>(std::ceil((this->getInteractionLength()) * _towerSideLengthReciprocal))
   */
  int _interactionLengthInTowers;

  /**
   * Specifies if the saved neighbors use newton 3 or not.
   */
  bool _neighborListIsNewton3;
};

}  // namespace autopas
