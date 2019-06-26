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
class VerletClusterLists : public ParticleContainer<Particle, FullParticleCell<Particle>> {
  /**
   * the index type to access the particle cells
   */
  using index_t = typename VerletClusterMaths::index_t;

 public:
  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly the
   * same side length. The rebuildFrequency should be chosen, s.t. the particles do
   * not move more than a distance of skin/2 between two rebuilds of the lists.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param skin the skin radius
   * @param rebuildFrequency specifies after how many pair-wise traversals the
   * neighbor lists are to be rebuild. A frequency of 1 means that they are
   * always rebuild, 10 means they are rebuild after 10 traversals.
   * @param clusterSize size of clusters
   */
  VerletClusterLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff,
                     double skin = 0, unsigned int rebuildFrequency = 1, int clusterSize = 4)
      : ParticleContainer<Particle, FullParticleCell<Particle>>(boxMin, boxMax, cutoff + skin),
        _clusterSize(clusterSize),
        _numClusters(0),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _skin(skin),
        _cutoff(cutoff),
        _cutoffSqr(cutoff * cutoff),
        _traversalsSinceLastRebuild(UINT_MAX),
        _rebuildFrequency(rebuildFrequency),
        _neighborListIsValid(false),
        _aosToSoaMapValid(false) {
    rebuild();
  }

  ContainerOption getContainerType() override { return ContainerOption::verletClusterLists; }

  /**
   * Function to iterate over all pairs of particles. (Only AoS)
   * This function only handles short-range interactions.
   * @tparam The type of the ParticleFunctor.
   * @tparam Traversal The type of the traversal.
   * @param f not used
   * @param traversal The traversal to use for the iteration
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal) {
    if (traversal->getUseNewton3()) {
      /// @todo implement newton3 for VerletClusterLists
      AutoPasLog(error, "Newton3 not implemented yet.");
      autopas::utils::ExceptionHandler::exception("VerletClusterLists does not support newton3 yet.");
    }

    if (needsRebuild()) {
      this->rebuild();
    }

    if (traversal->getDataLayout() == DataLayoutOption::soa && not _aosToSoaMapValid) {
      buildAosToSoaMap();
    }

    auto *traversalInterface = dynamic_cast<VerletClustersTraversalInterface<Particle> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setTraversalInfo(this);
      traversalInterface->initClusterTraversal();
      traversalInterface->traverseParticlePairs();
      traversalInterface->endClusterTraversal();
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletClusterLists::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }

    // we iterated, so increase traversal counter
    _traversalsSinceLastRebuild++;
  }

  /**
   * @copydoc VerletLists::addParticle()
   */
  void addParticle(Particle &p) override {
    _neighborListIsValid = false;
    // add particle somewhere, because lists will be rebuild anyways
    _clusters[0].addParticle(p);
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticle(Particle &haloParticle) override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.addHaloParticle not yet implemented.");
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.deleteHaloParticles not yet implemented.");
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  void updateContainer() override {
    AutoPasLog(debug, "updating container");
    _neighborListIsValid = false;
  }

  bool isContainerUpdateNeeded() override {
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.isContainerUpdateNeeded not yet implemented");
    return false;
  }

  TraversalSelectorInfo<FullParticleCell<Particle>> getTraversalSelectorInfo() override {
    return TraversalSelectorInfo<FullParticleCell<Particle>>(_cellsPerDim);
  }

  /**
   * Specifies whether the neighbor lists need to be rebuild.
   * @return true if the neighbor lists need to be rebuild, false otherwise
   */
  bool needsRebuild() {
    AutoPasLog(debug, "VerletLists: neighborlist is valid: {}", _neighborListIsValid);
    // if the neighbor list is NOT valid or we have not rebuild for _rebuildFrequency steps
    return (not _neighborListIsValid) or (_traversalsSinceLastRebuild >= _rebuildFrequency);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, FullParticleCell<Particle>>(&this->_clusters));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                      const std::array<double, 3> &higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    // @todo implement this if bounding boxes are here
    autopas::utils::ExceptionHandler::exception("VerletClusterLists.getRegionIterator not yet implemented.");
    return ParticleIteratorWrapper<Particle>();
  }

  /**
   * Helper method to iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @tparam inParallel If the iteration should be executed in parallel or sequential.
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
   * Returns the AoSToSoAMap for usage in the traversals of this container.
   * @return the AoSToSoAMap.
   */
  const auto &getAosToSoaMap() const { return _aosToSoaMap; }

  /**
   * Returns the number of clusters in this container.
   * @return The number of clusters in this container.
   */
  auto getNumClusters() const { return _numClusters; }

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
        auto &grid = _clusters[index];
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
        auto &grid = _clusters[index];
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
   * Recalculate grids and clusters,
   * build verlet lists and
   * pad clusters.
   */
  void rebuild() {
    _aosToSoaMapValid = false;

    std::vector<Particle> invalidParticles = collectParticlesAndClearClusters();

    auto boxSize = ArrayMath::sub(_boxMax, _boxMin);

    _gridSideLength = guessOptimalGridSideLength(invalidParticles.size(), boxSize);
    _gridSideLengthReciprocal = 1 / _gridSideLength;

    _cellsPerDim = calculateCellsPerDim(boxSize);
    // _cellsPerDim[2] is always 1
    index_t numCells = _cellsPerDim[0] * _cellsPerDim[1];

    // resize to number of grids
    _clusters.resize(numCells);
    _neighborLists.resize(numCells);

    sortParticlesIntoClusters(invalidParticles);

    // sort by last dimension and reserve space for dummy particles
    for (auto &cluster : _clusters) {
      cluster.sortByDim(2);

      size_t size = cluster.numParticles();
      size_t rest = size % _clusterSize;
      if (rest > 0) cluster.reserve(size + (_clusterSize - rest));
    }

    clearNeighborLists();

    updateVerletLists();
    // fill last cluster with dummy particles, such that each cluster is a multiple of _clusterSize
    padClusters();

    _numClusters = calculateNumClusters();
  }

  /**
   * Takes all particles from all clusters and returns them. Clusters are cleared.
   * @return All particles in the container.
   */
  std::vector<Particle> collectParticlesAndClearClusters() {
    std::vector<Particle> invalidParticles;
    for (auto &cluster : _clusters) {
      for (auto it = cluster.begin(); it.isValid(); ++it) {
        invalidParticles.push_back(*it);
      }
      cluster.clear();
    }
    return invalidParticles;
  }

  /**
   * Guesses the optimal grid side length.
   * @param numParticles The number of particles in the container.
   * @param boxSize The size of the domain.
   * @return an estimated optimal grid side length.
   */
  virtual double guessOptimalGridSideLength(size_t numParticles, std::array<double, 3> boxSize) const {
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
        _clusters[index].addParticle(particle);
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
   */
  void updateVerletLists() {
    const int boxRange = static_cast<int>(std::ceil((_cutoff + _skin) * _gridSideLengthReciprocal));

    const int gridMaxX = _cellsPerDim[0] - 1;
    const int gridMaxY = _cellsPerDim[1] - 1;
    // for all grids
    for (int yi = 0; yi <= gridMaxY; yi++) {
      for (int xi = 0; xi <= gridMaxX; xi++) {
        auto &iGrid = _clusters[VerletClusterMaths::index1D(xi, yi, _cellsPerDim)];
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

    // the neighbor list is now valid
    _neighborListIsValid = true;
    _traversalsSinceLastRebuild = 0;
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
        if (distXYsqr <= _cutoffSqr) {
          auto &jGrid = _clusters[VerletClusterMaths::index1D(xj, yj, _cellsPerDim)];
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
    float iBBoxBot = iGrid[iClusterIndex * _clusterSize].getR()[2];
    float iBBoxTop = iGrid[iClusterIndex * _clusterSize + iClusterSize - 1].getR()[2];
    auto &iClusterVerlet = iNeighbors[iClusterIndex];

    // iterate over full clusters of j-th grid.
    for (index_t zj = 0; zj < jSize; zj++) {
      addJClusterAsNeighborIfInRange(jGrid, zj, _clusterSize, iClusterVerlet, distXYsqr, iBBoxBot, iBBoxTop);
    }
    // special case: last cluster not full
    if (jRest > 0) {
      addJClusterAsNeighborIfInRange(jGrid, jSize, jRest, iClusterVerlet, distXYsqr, iBBoxBot, iBBoxTop);
    }
  }

  /**
   * Adds the given cluster in jGrid to the given neighbor list (iClusterNeighborList), if it is within the cutoff
   * radius.
   * @param jGrid The j-th grid.
   * @param jClusterIndex The index of the cluster to work on in the j-th grid.
   * @param jClusterSize The size of the cluster to work on in the j-th grid.
   * @param iClusterNeighborList The neighbor list of the cluster in the i-th grid to fill the neighbors for.
   * @param distXYsqr The distance between the i-th grid and the j-th grid in the xy-plane.
   * @param iBBoxBot The bottom z-coordinate of the cluster in the i-th grid.
   * @param iBBoxTop The top z-coordinate of the cluster in the i-th grid.
   */
  void addJClusterAsNeighborIfInRange(FullParticleCell<Particle> &jGrid, index_t jClusterIndex, int jClusterSize,
                                      std::vector<Particle *> &iClusterNeighborList, double distXYsqr, float iBBoxBot,
                                      float iBBoxTop) {
    Particle *jClusterStart = &jGrid[jClusterIndex * _clusterSize];
    // bbox in z of jGrid
    float jBBoxBot = jClusterStart->getR()[2];
    float jBBoxTop = (jClusterStart + (jClusterSize - 1))->getR()[2];

    double distZ = bboxDistance(iBBoxBot, iBBoxTop, jBBoxBot, jBBoxTop);
    if (distXYsqr + distZ * distZ <= _cutoffSqr) {
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
        auto &grid = _clusters[VerletClusterMaths::index1D(x, y, _cellsPerDim)];
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
   * Calculates the number of clusters in this container.
   * @return The number of clusters in this container.
   */
  index_t calculateNumClusters() {
    // iterate over all clusters
    index_t currentClusterIndex = 0;
    for (index_t x = 0; x < _cellsPerDim[0]; x++) {
      for (index_t y = 0; y < _cellsPerDim[1]; y++) {
        index_t index = VerletClusterMaths::index1D(x, y, _cellsPerDim);
        auto &grid = _clusters[index];

        const index_t numClustersInGrid = grid.numParticles() / _clusterSize;
        for (index_t clusterIndex = 0; clusterIndex < numClustersInGrid; clusterIndex++) {
          currentClusterIndex++;
        }
      }
    }
    return currentClusterIndex;
  }

  /**
   * Calculates the distance of two bounding boxes in one dimension.
   * @param min1 minimum coordinate of first bbox in tested dimension
   * @param max1 maximum coordinate of first bbox in tested dimension
   * @param min2 minimum coordinate of second bbox in tested dimension
   * @param max2 maximum coordinate of second bbox in tested dimension
   * @return distance
   */
  inline float bboxDistance(const float min1, const float max1, const float min2, const float max2) const {
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
   * Builds the _aosToSoaMap to be up to date with _clusters.
   */
  void buildAosToSoaMap() {
    index_t currentMapIndex = 0;

    const auto _clusterTraverseFunctor = [this, &currentMapIndex](Particle *clusterStart, int clusterSize,
                                                              std::vector<Particle *> &clusterNeighborList) {
      _aosToSoaMap[clusterStart] = currentMapIndex++;
    };

    // Cannot be parallelized at the moment.
    this->template traverseClusters<false>(_clusterTraverseFunctor);

    _aosToSoaMapValid = true;
  }

 private:
  /// neighbors of clusters for each grid
  std::vector<std::vector<std::vector<Particle *>>> _neighborLists;

  /// internal storage, particles are split into a grid in xy-dimension
  std::vector<FullParticleCell<Particle>> _clusters;
  int _clusterSize;

  /// The number of clusters. This is not equal to _cluseters.size(), as every grid might contain multiple clusters.
  index_t _numClusters;

  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;

  /// side length of xy-grid and reciprocal
  double _gridSideLength;
  double _gridSideLengthReciprocal;

  /// dimensions of grid
  std::array<index_t, 3> _cellsPerDim;

  /// skin radius
  double _skin;

  /// cutoff
  double _cutoff;
  double _cutoffSqr;

  /// how many pairwise traversals have been done since the last traversal
  unsigned int _traversalsSinceLastRebuild;

  /// specifies after how many pairwise traversals the neighbor list is to be
  /// rebuild
  unsigned int _rebuildFrequency;

  /// specifies if the neighbor list is currently valid
  bool _neighborListIsValid;

  /// Maps indices to the starting pointers for each cluster
  std::unordered_map<Particle *, index_t> _aosToSoaMap;

  /// If _aosToSoaMap is valid currently; that means up to date with clusters.
  bool _aosToSoaMapValid;
};

}  // namespace autopas
