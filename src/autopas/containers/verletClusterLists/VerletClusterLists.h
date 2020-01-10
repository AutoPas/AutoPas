/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <cmath>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/verletClusterLists/VerletClusterMaths.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/utils/ArrayMath.h"
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
class VerletClusterLists : public ParticleContainer<FullParticleCell<Particle>> {
  /**
   * the index type to access the particle cells
   */
  typedef VerletClusterMaths::index_t index_t;

 public:
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
                     double skin = 0, int clusterSize = 4)
      : ParticleContainer<FullParticleCell<Particle>>(boxMin, boxMax, cutoff, skin),
        _clusterSize(clusterSize),
        _numClusters(0),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _skin(skin),
        _cutoff(cutoff),
        _neighborListIsNewton3(false),
        _interactionLengthSqr((cutoff + skin) * (cutoff + skin)) {
    rebuild(false);
  }

  ContainerOption getContainerType() const override { return ContainerOption::verletClusterLists; }

  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", traversal->getTraversalType().to_string());

    auto *traversalInterface = dynamic_cast<VerletClustersTraversalInterface<Particle> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setClusterLists(*this);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletClusterLists::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle &p) override {
    // add particle somewhere, because lists will be rebuild anyways
    this->_cells[0].addParticle(p);
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(Particle &haloParticle) override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.addHaloParticle not yet implemented.");
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(Particle &haloParticle) override { throw std::runtime_error("not yet implemented"); }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    // quick and dirty: iterate over all particles and delete halo particles
    // @todo: make this proper
    for (auto iter = this->begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
      if (not iter->isOwned()) {
        internal::deleteParticle(iter);
      }
    }
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<Particle> updateContainer() override {
    AutoPasLog(debug, "updating container");
    // first delete all particles
    this->deleteHaloParticles();

    // next find invalid particles
    std::vector<Particle> invalidParticles;
    /// @todo: parallelize
    for (auto iter = this->begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      if (not utils::inBox(iter->getR(), _boxMin, _boxMax)) {
        invalidParticles.push_back(*iter);
        internal::deleteParticle(iter);
      }
    }

    return invalidParticles;
  }

  bool isContainerUpdateNeeded() const override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.isContainerUpdateNeeded not yet implemented");
    return false;
  }

  TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(_cellsPerDim, this->getInteractionLength(), {0., 0., 0.}, _clusterSize);
  }

  ParticleIteratorWrapper<Particle, true> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle, true>(
        new internal::ParticleIterator<Particle, FullParticleCell<Particle>, true>(&this->_cells));
  }

  ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    return ParticleIteratorWrapper<Particle, false>(
        new internal::ParticleIterator<Particle, FullParticleCell<Particle>, false>(&this->_cells));
  }

  ParticleIteratorWrapper<Particle, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle, true>();
  }

  ParticleIteratorWrapper<Particle, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    // @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle, false>();
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override { rebuild(traversal->getUseNewton3()); }

  /**
   * Helper method to iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @tparam inParallel If the iteration should be executed in parallel or sequential.  See traverseClustersParallel()
   * for thread safety.
   * @param loopBody The lambda to execute for all clusters. Parameters given are Particle* clusterStart, int
   * clusterSize, std::vector<Particle*> clusterNeighborList.
   */
  template <bool inParallel, class LoopBody>
  void traverseClusters(LoopBody &&loopBody) {
    if (inParallel) {
      traverseClustersParallel<LoopBody>(std::forward<LoopBody>(loopBody));
    } else {
      traverseClustersSequential<LoopBody>(std::forward<LoopBody>(loopBody));
    }
  }

  /**
   * Returns the ClusterIndexMap for usage in the traversals of this container.
   * @return the ClusterIndexMap.
   */
  const auto &getClusterIndexMap() const { return _clusterIndexMap; }

  /**
   * Returns the number of clusters in this container.
   * @return The number of clusters in this container.
   */
  auto getNumClusters() const { return _numClusters; }

  /**
   * Returns the neighbor lists of this container.
   * @return the neighbor lists of this container.
   */
  const auto &getNeighborLists() const { return _neighborLists; }

  /**
   * Returns the grid side length of the grids in the container.
   * @return the grid side length of the grids in the container.
   */
  auto getGridSideLength() const { return _gridSideLength; }

  /**
   * Returns the number of grids per dimension on the container.
   * @return the number of grids per dimension on the container.
   */
  auto getCellsPerDimension() const { return _cellsPerDim; }

  /**
   * Returns the 2D grid for the XY-plane of this container that defines the cluster towers.
   * @return the grids of this container for usage in traversals.
   */
  auto &getGrids() { return this->_cells; }

  /**
   * Returns the number of particles in each cluster.
   * @return the number of particles in each cluster.
   */
  auto getClusterSize() const { return _clusterSize; }

 protected:
  /**
   * Helper method to sequentially iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given are Particle* clusterStart, index_t
   * clusterSize, std::vector<Particle*> clusterNeighborList.
   */
  template <class LoopBody>
  void traverseClustersSequential(LoopBody &&loopBody) {
    for (index_t x = 0; x < _cellsPerDim[0]; x++) {
      for (index_t y = 0; y < _cellsPerDim[1]; y++) {
        index_t index = VerletClusterMaths::index1D(x, y, _cellsPerDim);
        auto &grid = this->_cells[index];
        auto &gridNeighborList = _neighborLists[index];

        const index_t numClustersInGrid = grid.numParticles() / _clusterSize;
        for (index_t clusterInGrid = 0; clusterInGrid < numClustersInGrid; clusterInGrid++) {
          Particle *iClusterStart = &grid[clusterInGrid * _clusterSize];
          auto &clusterNeighborList = gridNeighborList[clusterInGrid];
          loopBody(iClusterStart, _clusterSize, clusterNeighborList);
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
   * @param loopBody The lambda to execute for all clusters. Parameters given are Particle* clusterStart, index_t
   * clusterSize, std::vector<Particle*> clusterNeighborList.
   */
  template <class LoopBody>
  void traverseClustersParallel(LoopBody &&loopBody) {
    const index_t endX = _cellsPerDim[0];
    const index_t endY = _cellsPerDim[1];
#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (index_t x = 0; x < endX; x++) {
      for (index_t y = 0; y < endY; y++) {
        index_t index = VerletClusterMaths::index1D(x, y, _cellsPerDim);
        auto &grid = this->_cells[index];
        auto &gridNeighborList = _neighborLists[index];

        const index_t numClustersInGrid = grid.numParticles() / _clusterSize;
        for (index_t clusterInGrid = 0; clusterInGrid < numClustersInGrid; clusterInGrid++) {
          Particle *iClusterStart = &grid[clusterInGrid * _clusterSize];
          auto &clusterNeighborList = gridNeighborList[clusterInGrid];
          loopBody(iClusterStart, _clusterSize, clusterNeighborList);
        }
      }
    }
  }

  /**
   * Recalculate grids and clusters, build verlet lists and pad clusters.
   * @param useNewton3 If the everything should be build using newton 3 or not.
   */
  void rebuild(bool useNewton3) {
    std::vector<Particle> invalidParticles = collectParticlesAndClearClusters();

    auto boxSize = utils::ArrayMath::sub(_boxMax, _boxMin);

    _gridSideLength = estimateOptimalGridSideLength(invalidParticles.size(), boxSize);
    _gridSideLengthReciprocal = 1 / _gridSideLength;

    _cellsPerDim = calculateCellsPerDim(boxSize);
    // _cellsPerDim[2] is always 1
    index_t numCells = _cellsPerDim[0] * _cellsPerDim[1];

    // resize to number of grids
    this->_cells.resize(numCells);
    _neighborLists.resize(numCells);

    sortParticlesIntoClusters(invalidParticles);

    // sort by last dimension and reserve space for dummy particles
    for (auto &cluster : this->_cells) {
      cluster.sortByDim(2);

      size_t size = cluster.numParticles();
      size_t rest = size % _clusterSize;
      if (rest > 0) cluster.reserve(size + (_clusterSize - rest));
    }

    clearNeighborLists();

    _numClusters = buildClusterIndexMap();

    updateVerletLists(useNewton3);
    // fill last cluster with dummy particles, such that each cluster is a multiple of _clusterSize
    padClusters();
  }

  /**
   * Takes all particles from all clusters and returns them. Clusters are cleared.
   * @return All particles in the container.
   */
  std::vector<Particle> collectParticlesAndClearClusters() {
    std::vector<Particle> invalidParticles;
    for (auto &cluster : this->_cells) {
      for (auto it = cluster.begin(); it.isValid(); ++it) {
        invalidParticles.push_back(*it);
      }
      cluster.clear();
    }
    return invalidParticles;
  }

  /**
   * Estimates the optimal grid side length.
   * @param numParticles The number of particles in the container.
   * @param boxSize The size of the domain.
   * @return an estimated optimal grid side length.
   */
  virtual double estimateOptimalGridSideLength(size_t numParticles, std::array<double, 3> boxSize) const {
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
   * Calculates the cells per dimension in the container using the _gridSideLengthReciprocal.
   * @param boxSize the size of the domain.
   * @return the cells per dimension in the container.
   */
  std::array<index_t, 3> calculateCellsPerDim(std::array<double, 3> boxSize) const {
    std::array<index_t, 3> cellsPerDim{};
    for (int d = 0; d < 2; d++) {
      cellsPerDim[d] = static_cast<index_t>(std::ceil(boxSize[d] * _gridSideLengthReciprocal));
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
  void sortParticlesIntoClusters(std::vector<Particle> &particles) {
    for (auto &particle : particles) {
      if (utils::inBox(particle.getR(), _boxMin, _boxMax)) {
        auto index = get1DIndexOfPosition(particle.getR());
        this->_cells[index].addParticle(particle);
      }
    }
  }

  /**
   * Clears all neighbor lists.
   */
  void clearNeighborLists() {
    for (auto &verlet : _neighborLists) {
      verlet.clear();
    }
  }

  /**
   * Update the verlet lists.
   *
   * @param useNewton3 If newton 3 should be used to build the neighbor lists or not. If true, only saves neighbor
   * clusters that have a higher index that the current cluster. (@see buildClusterIndexMap())
   */
  void updateVerletLists(bool useNewton3) {
    _neighborListIsNewton3 = useNewton3;

    const int boxRange = static_cast<int>(std::ceil((_cutoff + _skin) * _gridSideLengthReciprocal));

    const int gridMaxX = _cellsPerDim[0] - 1;
    const int gridMaxY = _cellsPerDim[1] - 1;
    // for all grids
    for (int yi = 0; yi <= gridMaxY; yi++) {
      for (int xi = 0; xi <= gridMaxX; xi++) {
        auto &iGrid = this->_cells[VerletClusterMaths::index1D(xi, yi, _cellsPerDim)];
        // calculate number of full clusters and rest
        index_t iSize = iGrid.numParticles() / _clusterSize;
        int iRest = iGrid.numParticles() % _clusterSize;

        const int minX = std::max(xi - boxRange, 0);
        const int minY = std::max(yi - boxRange, 0);
        const int maxX = std::min(xi + boxRange, gridMaxX);
        const int maxY = std::min(yi + boxRange, gridMaxY);

        auto &iNeighbors = _neighborLists[VerletClusterMaths::index1D(xi, yi, _cellsPerDim)];
        if (iRest > 0)
          iNeighbors.resize(iSize + 1);
        else
          iNeighbors.resize(iSize);

        addClustersOfNeighborGridsAsNeighborsIfInRange(iGrid, iSize, iRest, iNeighbors, minX, maxX, minY, maxY, xi, yi);
      }
    }
  }

  /**
   * Iterates over neighbor grids of the i-th grid and adds all clusters in them that are within the cutoff radius to
   * the neighbor list of the clusters in the i-th grid.
   * @param iGrid The i-th grid.
   * @param iSize The number of full clusters in the i-th grid.
   * @param iRest If the last cluster is not full: The number of particles in the last cluster. 0 otherwise.
   * @param iNeighbors The neighbor list of the i-th grid.
   * @param minX
   * @param maxX
   * @param minY
   * @param maxY
   * @param xi The x-index of the i-th grid.
   * @param yi the y-index of the i-th grid.
   */
  void addClustersOfNeighborGridsAsNeighborsIfInRange(FullParticleCell<Particle> &iGrid, index_t iSize, int iRest,
                                                      std::vector<std::vector<Particle *>> &iNeighbors, const int minX,
                                                      const int maxX, const int minY, const int maxY, const int xi,
                                                      const int yi) {
    // for all neighbor grids
    for (int yj = minY; yj <= maxY; yj++) {
      double distY = std::max(0, std::abs(yi - yj) - 1) * _gridSideLength;

      for (int xj = minX; xj <= maxX; xj++) {
        double distX = std::max(0, std::abs(xi - xj) - 1) * _gridSideLength;

        // calculate distance in xy-plane and skip if already longer than cutoff
        double distXYsqr = distX * distX + distY * distY;
        if (distXYsqr <= _interactionLengthSqr) {
          auto &jGrid = this->_cells[VerletClusterMaths::index1D(xj, yj, _cellsPerDim)];
          // calculate number of  full clusters and rest
          const index_t jSize = jGrid.numParticles() / _clusterSize;
          const int jRest = jGrid.numParticles() % _clusterSize;

          // for all clusters in the i-th grid
          for (index_t zi = 0; zi < iSize; zi++) {
            addAllJClustersAsNeighborIfInRange(iGrid, zi, _clusterSize, iNeighbors, jGrid, jSize, jRest, distXYsqr);
          }

          // special case: last cluster of iGrid not full
          if (iRest > 0) {
            addAllJClustersAsNeighborIfInRange(iGrid, iSize, iRest, iNeighbors, jGrid, jSize, jRest, distXYsqr);
          }
        }
      }
    }
  }

  /**
   * Adds all clusters in jGrid that are within the cutoff radius to the neighbor list of the given cluster in iGrid
   * (iClusterIndex).
   * @param iGrid The i-th grid.
   * @param iClusterIndex The index of the cluster to work on in the i-th grid.
   * @param iClusterSize The size of th cluster with index iClusterIndex in the i-th grid.
   * @param iNeighbors The neighbor list of the i-th grid.
   * @param jGrid The j-th grid.
   * @param jSize The number of full clusters in the j-th grid.
   * @param jRest If the last cluster is not full: The number of particles in the last cluster. 0 otherwise.
   * @param distXYsqr The distance between the i-th grid and the j-th grid in the xy-plane.
   */
  void addAllJClustersAsNeighborIfInRange(FullParticleCell<Particle> &iGrid, index_t iClusterIndex, int iClusterSize,
                                          std::vector<std::vector<Particle *>> &iNeighbors,
                                          FullParticleCell<Particle> &jGrid, index_t jSize, int jRest,
                                          double distXYsqr) {
    // bbox in z of iGrid
    double iBBoxBot = iGrid[iClusterIndex * _clusterSize].getR()[2];
    double iBBoxTop = iGrid[iClusterIndex * _clusterSize + iClusterSize - 1].getR()[2];
    auto &iClusterNeighborList = iNeighbors[iClusterIndex];
    Particle *iClusterStart = &iGrid[iClusterIndex * _clusterSize];

    // iterate over full clusters of j-th grid.
    for (index_t jClusterIndex = 0; jClusterIndex < jSize; jClusterIndex++) {
      Particle *jClusterStart = &jGrid[jClusterIndex * _clusterSize];
      // If newton 3 is used, only add clusters as neighbors that have a equal or higher index. Skip otherwise.
      if (_neighborListIsNewton3 and _clusterIndexMap.at(iClusterStart) > _clusterIndexMap.at(jClusterStart)) continue;

      addJClusterAsNeighborIfInRange(jGrid, jClusterStart, _clusterSize, iClusterNeighborList, distXYsqr, iBBoxBot,
                                     iBBoxTop);
    }
    // special case: last cluster not full
    if (jRest > 0) {
      Particle *jClusterStart = &jGrid[jSize * _clusterSize];
      // If newton 3 is used, only add clusters as neighbors that have a equal or higher index. Skip otherwise.
      if (not(_neighborListIsNewton3 and _clusterIndexMap.at(iClusterStart) > _clusterIndexMap.at(jClusterStart))) {
        addJClusterAsNeighborIfInRange(jGrid, jClusterStart, jRest, iClusterNeighborList, distXYsqr, iBBoxBot,
                                       iBBoxTop);
      }
    }
  }

  /**
   * Adds the given cluster in jGrid to the given neighbor list (iClusterNeighborList), if it is within the cutoff
   * radius.
   * @param jGrid The j-th grid.
   * @param jClusterStart A pointer to the start of the cluster to work on in the j-th grid.
   * @param jClusterSize The size of the cluster to work on in the j-th grid.
   * @param iClusterNeighborList The neighbor list of the cluster in the i-th grid to fill the neighbors for.
   * @param distXYsqr The distance between the i-th grid and the j-th grid in the xy-plane.
   * @param iBBoxBot The bottom z-coordinate of the cluster in the i-th grid.
   * @param iBBoxTop The top z-coordinate of the cluster in the i-th grid.
   */
  void addJClusterAsNeighborIfInRange(FullParticleCell<Particle> &jGrid, Particle *jClusterStart, int jClusterSize,
                                      std::vector<Particle *> &iClusterNeighborList, double distXYsqr, double iBBoxBot,
                                      double iBBoxTop) {
    // bbox in z of jGrid
    double jBBoxBot = jClusterStart->getR()[2];
    double jBBoxTop = (jClusterStart + (jClusterSize - 1))->getR()[2];

    double distZ = bboxDistance(iBBoxBot, iBBoxTop, jBBoxBot, jBBoxTop);
    if (distXYsqr + distZ * distZ <= _interactionLengthSqr) {
      iClusterNeighborList.push_back(jClusterStart);
    }
  }

  /**
   * Pad clusters with dummy particles
   * until each cluster is a multiple of _clusterSize.
   * Useful for SIMD vectorization.
   */
  void padClusters() {
    for (index_t x = 0; x < _cellsPerDim[0]; x++) {
      for (index_t y = 0; y < _cellsPerDim[1]; y++) {
        auto &grid = this->_cells[VerletClusterMaths::index1D(x, y, _cellsPerDim)];
        index_t rest = grid.numParticles() % _clusterSize;
        if (rest > 0) {
          for (int i = rest; i < _clusterSize; i++) {
            Particle p = Particle();
            p.setR({2 * x * _cutoff, 2 * y * _cutoff, 2 * _boxMax[2] + 2 * i * _cutoff});
            grid.addParticle(p);
          }
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
  inline double bboxDistance(const double min1, const double max1, const double min2, const double max2) const {
    if (max1 < min2) {
      return min2 - max1;
    } else if (min1 > max2) {
      return min1 - max2;
    } else {
      return 0;
    }
  }

  /**
   * Gets the 1d grid index containing a particle in given position.
   * @param pos the position of the particle
   * @return the index of the grid
   */
  inline index_t get1DIndexOfPosition(const std::array<double, 3> &pos) const {
    std::array<index_t, 2> cellIndex{};

    for (int dim = 0; dim < 2; dim++) {
      const long int value = (static_cast<long int>(floor((pos[dim] - _boxMin[dim]) * _gridSideLengthReciprocal))) + 1l;
      const index_t nonnegativeValue = static_cast<index_t>(std::max(value, 0l));
      const index_t nonLargerValue = std::min(nonnegativeValue, _cellsPerDim[dim] - 1);
      cellIndex[dim] = nonLargerValue;
      /// @todo this is a sanity check to prevent doubling of particles, but
      /// could be done better! e.g. by border and flag manager
      if (pos[dim] >= _boxMax[dim]) {
        cellIndex[dim] = _cellsPerDim[dim] - 1;
      } else if (pos[dim] < _boxMin[dim]) {
        cellIndex[dim] = 0;
      }
    }

    return VerletClusterMaths::index1D(cellIndex[0], cellIndex[1], _cellsPerDim);
  }

  /**
   * Builds the _clusterIndexMap to be up to date with _cells.
   *
   * Every cluster gets an index assigned. The indices are given in a way so that the VerletClustersColoringTraversal
   * works as easy as possible with newton 3. The newton 3 neighbor list just has to only save neighbors with a higher
   * index, and there will be no data races.
   *
   * For each cluster now holds (with x-axis as left <=> right, y-axis <=> as top <=> bottom):
   *    - The indices of all clusters of the three color cells above and the color cell to the left are lower.
   *    - The indices of all clusters of the three color cells below and the color cell to the right are higher.
   *    - For all grids of the same color cell holds:
   *        - The indices of all clusters of the three grids above and the grids to the left are lower.
   *        - The indices of all clusters of the three grids below and the grids to the right are higher.
   *    - For all clusters in the same grid holds:
   *        - The indices of all clusters with a lower z-coordinate than the current cluster are lower.
   *        - The indices of all clusters with a higher z-coordinate than the current cluster are higher.
   *
   * @return The number of clusters in the container.
   */
  index_t buildClusterIndexMap() {
    index_t nextFreeMapIndex = 0;

    int gridsPerColoringCell = static_cast<int>(std::ceil((_cutoff + _skin) / _gridSideLength));
    std::array<unsigned long, 3> coloringCellsPerDim{};
    for (int i = 0; i < 3; i++) {
      coloringCellsPerDim[i] =
          static_cast<unsigned long>(std::ceil(_cellsPerDim[i] / static_cast<double>(gridsPerColoringCell)));
    }

    for (unsigned long yColorCell = 0; yColorCell < coloringCellsPerDim[1]; yColorCell++) {
      for (unsigned long xColorCell = 0; xColorCell < coloringCellsPerDim[0]; xColorCell++) {
        nextFreeMapIndex = indexColorCell(xColorCell, yColorCell, gridsPerColoringCell, nextFreeMapIndex);
      }
    }

    return nextFreeMapIndex;
  }

 private:
  /**
   * Indexes all clusters of one color cell (inserts value into _clusterIndexMap) starting with currentMapIndex.
   *
   * The scheme follows the documentation from buildClusterIndexMap().
   * @param xColorCell The x coordinate of the color cell.
   * @param yColorCell The y coordinate of the color cell.
   * @param gridsPerColoringCell The number of grids in x and y dimension of this color cell.
   * @param currentMapIndex The first index to use.
   * @return The next available index after this cell.
   */
  index_t indexColorCell(unsigned long xColorCell, unsigned long yColorCell, int gridsPerColoringCell,
                         index_t currentMapIndex) {
    for (int yInner = 0; yInner < gridsPerColoringCell; yInner++) {
      for (int xInner = 0; xInner < gridsPerColoringCell; xInner++) {
        unsigned long y = yColorCell * gridsPerColoringCell + yInner;
        unsigned long x = xColorCell * gridsPerColoringCell + xInner;

        // Not every coloring cell has to have gridsPerColoringCell grids in every direction.
        if (x >= _cellsPerDim[0] or y >= _cellsPerDim[1]) {
          continue;
        }

        unsigned long gridIndex1D = VerletClusterMaths::index1D(x, y, _cellsPerDim);
        auto &currentGrid = this->_cells[gridIndex1D];

        auto numClusters = currentGrid.numParticles() / _clusterSize;
        int rest = currentGrid.numParticles() % _clusterSize;
        if (rest > 0) numClusters++;

        for (unsigned long currentCluster = 0; currentCluster < numClusters; currentCluster++) {
          Particle *clusterStart = &currentGrid[currentCluster * _clusterSize];
          _clusterIndexMap[clusterStart] = currentMapIndex++;
        }
      }
    }
    return currentMapIndex;
  }

 private:
  /**
   * Neighbors of clusters for each grid. If it uses newton 3 is saved in _neighborListIsNewton3.
   * If it uses newton 3: Only the neighbor clusters that have a higher index are saved. (@see _clusterIndexMap)
   */
  std::vector<std::vector<std::vector<Particle *>>> _neighborLists;

  /**
   * The number of particles in a full cluster.
   */
  int _clusterSize;

  /**
   * The number of clusters. This is not equal to _cells.size(), as every grid (=cell) might contain multiple clusters.
   */
  index_t _numClusters;

  /**
   * Box min of the domain.
   */
  std::array<double, 3> _boxMin;

  /**
   * Box max of the domain.
   */
  std::array<double, 3> _boxMax;

  /**
   * Side length of xy-grid.
   */
  double _gridSideLength{0.};

  /**
   *  Reciprocal of _gridSideLength.
   */
  double _gridSideLengthReciprocal{0.};

  /**
   * Dimensions of the 2D xy-grid.
   */
  std::array<index_t, 3> _cellsPerDim{};

  /**
   * The skin radius.
   */
  double _skin;

  /**
   * The cutoff.
   */
  double _cutoff;

  /**
   * Specifies if the neighbor list uses newton 3 or not.
   */
  bool _neighborListIsNewton3;

  /**
   * Maps indices to the starting pointers for each cluster. For the idea behind the assignment, @see
   * buildClusterIndexMap().
   */
  std::unordered_map<Particle *, index_t> _clusterIndexMap;

  /**
   * (_cutoff + _skin)^2.
   */
  double _interactionLengthSqr;
};

}  // namespace autopas
