/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <cmath>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleDeletedObserver.h"
#include "autopas/containers/UnknowingCellBorderAndFlagManager.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "autopas/containers/verletClusterLists/VerletClusterListsRebuilder.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Timer.h"

namespace autopas {

/**
 * Particles are divided into clusters.
 * The VerletClusterLists class uses neighborhood lists for each cluster
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * @tparam Particle
 */
template <class Particle>
class VerletClusterLists : public ParticleContainerInterface<FullParticleCell<Particle>>,
                           public internal::ParticleDeletedObserver {
 public:
  /**
   * Defines a cluster range used in the static cluster-thread-partition.
   */
  struct ClusterRange {
    /**
     * The index of the tower that contains the first cluster.
     */
    size_t startTowerIndex{};
    /**
     * The index of the first cluster in its tower.
     */
    size_t startIndexInTower{};
    /**
     * The number of clusters in the range.
     */
    size_t numClusters{};
  };

  /**
   * Constructor of the VerletClusterLists class.
   * The neighbor lists are build using a estimated density.
   * The box is divided into cuboids with roughly the
   * same side length.
   * @param boxMin The lower corner of the domain.
   * @param boxMax The upper corner of the domain.
   * @param cutoff The cutoff radius of the interaction.
   * @param skin The skin radius.
   * @param clusterSize Number of particles per cluster.
   */
  VerletClusterLists(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff, double skin,
                     size_t clusterSize)
      : ParticleContainerInterface<FullParticleCell<Particle>>(),
        _clusterSize{clusterSize},
        _numClusters{0},
        _numTowersPerInteractionLength{0},
        _boxMin{boxMin},
        _boxMax{boxMax},
        _haloBoxMin{utils::ArrayMath::subScalar(boxMin, cutoff + skin)},
        _haloBoxMax{utils::ArrayMath::addScalar(boxMax, cutoff + skin)},
        _cutoff{cutoff},
        _skin{skin} {
    // always have at least one tower.
    _towers.push_back(internal::ClusterTower<Particle>(_clusterSize));
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletClusterLists; }

  /**
   * @copydoc ParticleContainerInterface::iteratePairwise()
   */
  void iteratePairwise(TraversalInterface *traversal) override {
    if (_isValid == ValidityState::cellsAndListsValid) {
      autopas::utils::ExceptionHandler::exception(
          "VerletClusterLists::iteratePairwise(): Trying to do a pairwise iteration, even though verlet lists are not "
          "valid.");
    }
    auto *traversalInterface = dynamic_cast<VCLTraversalInterface<Particle> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setClusterLists(*this);
      traversalInterface->setTowers(_towers);
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
   * Adds the given particle to the container. rebuildVerletLists() has to be called to have it actually sorted in.
   * @param p The particle to add.
   */
  void addParticleImpl(const Particle &p) override {
    _isValid = ValidityState::invalid;
    _particlesToAdd.push_back(p);
  }

  /**
   * @copydoc VerletLists::addHaloParticle()
   */
  void addHaloParticleImpl(const Particle &haloParticle) override {
    _isValid = ValidityState::invalid;
    Particle copy = haloParticle;
    copy.setOwnershipState(OwnershipState::halo);
    _particlesToAdd.push_back(copy);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const Particle &haloParticle) override {
    Particle pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);

    for (auto it = getRegionIterator(utils::ArrayMath::subScalar(pCopy.getR(), this->getSkin() / 2),
                                     utils::ArrayMath::addScalar(pCopy.getR(), this->getSkin() / 2),
                                     IteratorBehavior::haloOnly);
         it.isValid(); ++it) {
      if (pCopy.getID() == it->getID()) {
        *it = pCopy;
        return true;
      }
    }
    return false;
  }

  /**
   * @copydoc VerletLists::deleteHaloParticles
   */
  void deleteHaloParticles() override {
    bool deletedSth = false;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(|| : deletedSth)
#endif
    {
      for (auto iter = this->begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
        internal::deleteParticle(iter);
        deletedSth = true;
      }
    }
    if (deletedSth) {
      _isValid = ValidityState::invalid;
    }
  }

  /**
   * @copydoc VerletLists::updateContainer()
   */
  [[nodiscard]] std::vector<Particle> updateContainer() override {
    // First delete all halo particles.
    this->deleteHaloParticles();
    // Delete dummy particles.
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
    for (size_t i = 0ul; i < _towers.size(); ++i) {
      _towers[i].deleteDummyParticles();
    }

    // next find invalid particles
    std::vector<Particle> invalidParticles;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
    {
      std::vector<Particle> myInvalidParticles;
      for (auto iter = this->begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
        if (not utils::inBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
          myInvalidParticles.push_back(*iter);
          internal::deleteParticle(iter);
        }
      }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
      invalidParticles.insert(invalidParticles.end(), myInvalidParticles.begin(), myInvalidParticles.end());
    }
    if (not invalidParticles.empty()) {
      _isValid = ValidityState::invalid;
    }
    return invalidParticles;
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    std::array<double, 3> towerSize = {_towerSideLength, _towerSideLength,
                                       this->getHaloBoxMax()[2] - this->getHaloBoxMin()[2]};
    std::array<unsigned long, 3> towerDimensions = {_towersPerDim[0], _towersPerDim[1], 1};
    return TraversalSelectorInfo(towerDimensions, this->getInteractionLength(), towerSize, _clusterSize);
  }

  /**
   * @copydoc ParticleContainerInterface::begin()
   * @note This function additionally rebuilds the towers if the tower-structure isn't valid.
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, true> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // For good openmp scalability we want the particles to be sorted into the clusters, so we do this!
#ifdef AUTOPAS_OPENMP
#pragma omp single
#endif
    if (_isValid == ValidityState::invalid) {
      rebuildTowersAndClusters();
    }
    // there is an implicit barrier at end of single!
    return ParticleIteratorWrapper<Particle, true>(
        new internal::ParticleIterator<Particle, internal::ClusterTower<Particle>, true>(
            &(this->_towers), 0, &unknowingCellBorderAndFlagManager, behavior));
  }

  /**
   * @copydoc ParticleContainerInterface::begin()
   * @note const version.
   * @note This function additionally iterates over the _particlesToAdd vector if the tower-structure isn't valid.
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    /// @todo use proper cellBorderAndFlagManager instead of the unknowing.
    if (_isValid != ValidityState::invalid) {
      if (not _particlesToAdd.empty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::begin() const: Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers.
      return ParticleIteratorWrapper<Particle, false>{
          new internal::ParticleIterator<Particle, internal::ClusterTower<Particle>, false>(
              &(this->_towers), 0, &unknowingCellBorderAndFlagManager, behavior)};
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      return ParticleIteratorWrapper<Particle, false>{
          new internal::ParticleIterator<Particle, internal::ClusterTower<Particle>, false>(
              &(this->_towers), 0, &unknowingCellBorderAndFlagManager, behavior, &_particlesToAdd)};
    }
  }

  /**
   * @copydoc ParticleContainerInterface::getRegionIterator()
   * @note This function additionally rebuilds the towers if the tower-structure isn't valid.
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // Special iterator requires sorted cells.
    // Only one thread is allowed to rebuild the towers, so we do an omp single here.
#ifdef AUTOPAS_OPENMP
#pragma omp single
#endif
    if (_isValid == ValidityState::invalid) {
      rebuildTowersAndClusters();
    }
    // there is an implicit barrier at end of single!

    auto [lowerCornerInBounds, upperCornerInBounds, cellsOfInterest] =
        getRegionIteratorHelper(lowerCorner, higherCorner, behavior);

    return ParticleIteratorWrapper<Particle, true>(
        new internal::RegionParticleIterator<Particle, internal::ClusterTower<Particle>, true>(
            &this->_towers, lowerCornerInBounds, upperCornerInBounds, cellsOfInterest,
            &internal::unknowingCellBorderAndFlagManager, behavior));
  }

  /**
   * @copydoc ParticleContainerInterface::getRegionIterator()
   * @note const version.
   * @note This function additionally iterates over _particlesToAdd if the container structure isn't valid.
   */
  [[nodiscard]] ParticleIteratorWrapper<Particle, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    if (_isValid != ValidityState::invalid && not _particlesToAdd.empty()) {
      autopas::utils::ExceptionHandler::exception(
          "VerletClusterLists::begin() const: Error: particle container is valid, but _particlesToAdd isn't empty!");
    }

    auto [lowerCornerInBounds, upperCornerInBounds, cellsOfInterest] =
        getRegionIteratorHelper(lowerCorner, higherCorner, behavior);

    return ParticleIteratorWrapper<Particle, false>(
        new internal::RegionParticleIterator<Particle, internal::ClusterTower<Particle>, false>(
            &this->_towers, lowerCornerInBounds, upperCornerInBounds, cellsOfInterest,
            &internal::unknowingCellBorderAndFlagManager, behavior,
            _isValid != ValidityState::invalid ? nullptr : &_particlesToAdd));
  }

  /**
   * @copydoc ParticleContainerInterface::rebuildNeighborLists()
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    if (_isValid == ValidityState::invalid) {
      rebuildTowersAndClusters();
    }
    _builder->rebuildNeighborListsAndFillClusters(traversal->getUseNewton3());

    auto *clusterTraversalInterface = dynamic_cast<VCLTraversalInterface<Particle> *>(traversal);
    if (clusterTraversalInterface) {
      if (clusterTraversalInterface->needsStaticClusterThreadPartition()) {
        calculateClusterThreadPartition();
      }
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletClusterLists::rebuildNeighborLists. TraversalID: {}",
          traversal->getTraversalType());
    }
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

  /**
   * @copydoc ParticleContainerInterface::getNumParticles()
   */
  [[nodiscard]] unsigned long getNumParticles() const override {
    unsigned long sum = 0;
    for (size_t index = 0; index < _towers.size(); index++) {
      sum += _towers[index].getNumActualParticles();
    }
    sum += _particlesToAdd.size();
    return sum;
  }

  /**
   * Returns the cluster-thread-partition.
   * @return The cluster-thread-partition.
   */
  const auto &getClusterThreadPartition() const { return _clusterThreadPartition; }

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
   * Returns 1 / towerSideLength
   * @return
   */
  auto getTowerSideLengthReciprocal() const { return _towerSideLengthReciprocal; }

  /**
   * Returns the number of grids per dimension on the container.
   * @return the number of grids per dimension on the container.
   */
  auto getTowersPerDimension() const { return _towersPerDim; }

  /**
   * Returns the number of particles in each cluster.
   * @return the number of particles in each cluster.
   */
  auto getClusterSize() const { return _clusterSize; }

  /**
   * Returns the towers per interaction length. That is how many towers fit into one interaction length rounded up.
   * @return the number of towers per interaction length.
   */
  auto getNumTowersPerInteractionLength() const { return _numTowersPerInteractionLength; }

  /**
   * Loads all particles of the container in their correct SoA and generates the SoAViews for the clusters.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for loading the particles into the SoA.
   */
  template <class Functor>
  void loadParticlesIntoSoAs(Functor *functor) {
    const auto numTowers = _towers.size();
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t index = 0; index < numTowers; index++) {
      _towers[index].loadSoA(functor);
    }
  }

  /**
   * Extracts all SoAs of the container into the particles.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for extracting the SoAs into the particles..
   */
  template <class Functor>
  void extractParticlesFromSoAs(Functor *functor) {
    const auto numTowers = _towers.size();
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t index = 0; index < numTowers; index++) {
      _towers[index].extractSoA(functor);
    }
  }

  /**
   * Returns a reference to the tower for the given tower grid coordinates.
   * @param x The x-th tower in x direction.
   * @param y The y-th tower in y direction.
   * @return a reference to the tower for the given tower grid coordinates.
   */
  auto &getTowerAtCoordinates(const size_t x, const size_t y) { return _towers[towerIndex2DTo1D(x, y)]; }

  /**
   * Returns the 1D index for the given tower grid coordinates of a tower.
   *
   * @param x The x-coordinate of the tower.
   * @param y The y-coordinate of the tower.
   * @param towersPerDim The number of towers in each dimension.
   * @return the 1D index for the given tower grid coordinates of a tower.
   */
  static auto towerIndex2DTo1D(const size_t x, const size_t y, const std::array<size_t, 2> towersPerDim) {
    return x + y * towersPerDim[0];
  }

  /**
   * Returns the 1D index for the given 2D-coordinates of a tower.
   *
   * @param x The x-coordinate of the tower.
   * @param y The y-coordinate of the tower.
   * @return the 1D index for the given 2D-coordinates of a tower.
   */
  [[nodiscard]] size_t towerIndex2DTo1D(const size_t x, const size_t y) const {
    return towerIndex2DTo1D(x, y, _towersPerDim);
  }

  /**
   * @copydoc ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override { return _boxMax; }

  /**
   * @copydoc ParticleContainerInterface::setBoxMax()
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override { _boxMax = boxMax; }

  /**
   * Get the upper corner of the halo box.
   * @return the upper corner of the halo box.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMax() const { return _haloBoxMax; }

  /**
   * @copydoc ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override { return _boxMin; }

  /**
   * @copydoc ParticleContainerInterface::setBoxMin()
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override { _boxMin = boxMin; }

  /**
   * Get the lower corner of the halo box.
   * @return the lower corner of the halo box.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMin() const { return _haloBoxMin; }

  /**
   * @copydoc ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const override { return _cutoff; }

  /**
   * @copydoc ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) override { _cutoff = cutoff; }

  /**
   * @copydoc ParticleContainerInterface::getSkin()
   */
  [[nodiscard]] double getSkin() const override { return _skin; }

  /**
   * @copydoc ParticleContainerInterface::setSkin()
   */
  void setSkin(double skin) override { _skin = skin; }

  /**
   * @copydoc ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const override { return _cutoff + _skin; }

  /**
   * @copydoc ParticleContainerInterface::deleteAllParticles()
   */
  void deleteAllParticles() override {
    _isValid = ValidityState::invalid;
    _particlesToAdd.clear();
    std::for_each(_towers.begin(), _towers.end(), [](auto &tower) { tower.clear(); });
  }

 protected:
  /**
   * Rebuild the towers and the clusters.
   * This function sets the container structure to valid.
   */
  void rebuildTowersAndClusters() {
    _builder = std::make_unique<internal::VerletClusterListsRebuilder<Particle>>(*this, _towers, _particlesToAdd,
                                                                                 _clusterSize);
    std::tie(_towerSideLength, _numTowersPerInteractionLength, _towersPerDim, _numClusters) =
        _builder->rebuildTowersAndClusters();
    _towerSideLengthReciprocal = 1 / _towerSideLength;
    _isValid = ValidityState::cellsValidListsInvalid;
    for (auto &tower : _towers) {
      tower.setParticleDeletionObserver(this);
    }
  }

  /**
   * Helper method to sequentially iterate over all clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <class LoopBody>
  void traverseClustersSequential(LoopBody &&loopBody) {
    for (size_t x = 0; x < _towersPerDim[0]; x++) {
      for (size_t y = 0; y < _towersPerDim[1]; y++) {
        auto &tower = getTowerAtCoordinates(x, y);
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
    const auto towersPerDimX = _towersPerDim[0];
    const auto towersPerDimY = _towersPerDim[1];
#if defined(AUTOPAS_OPENMP)
    /// @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (size_t x = 0; x < towersPerDimX; x++) {
      for (size_t y = 0; y < towersPerDimY; y++) {
        auto &tower = getTowerAtCoordinates(x, y);

        for (auto &cluster : tower.getClusters()) {
          loopBody(cluster);
        }
      }
    }
  }

  /**
   * Calculates a cluster thread partition that aims to give each thread about the same amount of cluster pair
   * interactions, if each thread handles the neighbors of all clusters it gets assigned.
   */
  void calculateClusterThreadPartition() {
    size_t numClusterPairs = 0;
    this->template traverseClusters<false>(
        [&numClusterPairs](auto &cluster) { numClusterPairs += cluster.getNeighbors().size(); });

    constexpr int minNumClusterPairsPerThread = 1000;
    auto numThreads =
        std::clamp(static_cast<int>(numClusterPairs / minNumClusterPairsPerThread), 1, autopas_get_max_threads());

    size_t numClusterPairsPerThread =
        std::max(static_cast<unsigned long>(std::ceil(static_cast<double>(numClusterPairs) / numThreads)), 1ul);
    if (numClusterPairsPerThread * numThreads < numClusterPairs) {
      autopas::utils::ExceptionHandler::exception(
          "VerletClusterLists::calculateClusterThreadPartition(): numClusterPairsPerThread ({}) * numThreads ({})={} "
          "should always "
          "be at least the amount of Cluster Pairs ({})!",
          numClusterPairsPerThread, numThreads, numClusterPairsPerThread * numThreads, numClusterPairs);
    }
    fillClusterRanges(numClusterPairsPerThread, numThreads);
  }

  /**
   * Fills in the cluster ranges of the cluster thread partition. It aims to assign each thread appropriately the same
   * number of cluster pairs.
   * @param numClusterPairsPerThread The approximate number of cluster pairs per thread.
   * @param numThreads The number of threads to use.
   */
  void fillClusterRanges(size_t numClusterPairsPerThread, int numThreads) {
    if (numClusterPairsPerThread < 1) {
      autopas::utils::ExceptionHandler::exception(
          "VerletClusterLists::fillClusterRanges(): numClusterPairsPerThread({}) is less than one, this is not "
          "supported "
          "and will lead to errors!",
          numClusterPairsPerThread);
    }
    _clusterThreadPartition.resize(numThreads);

    size_t currentThread = 0;
    size_t currentNumClustersToAdd = 0;
    size_t numClusterPairsTotal = 0;
    bool threadIsInitialized = false;
    // Iterate over the clusters of all towers
    for (size_t currentTowerIndex = 0; currentTowerIndex < _towers.size(); currentTowerIndex++) {
      auto &currentTower = _towers[currentTowerIndex];
      for (size_t currentClusterInTower = 0; currentClusterInTower < currentTower.getNumClusters();
           currentClusterInTower++) {
        auto &currentCluster = currentTower.getCluster(currentClusterInTower);

        // If on a new thread, start with the clusters for this thread here.
        if (not threadIsInitialized) {
          _clusterThreadPartition[currentThread] = {currentTowerIndex, currentClusterInTower, 0};
          threadIsInitialized = true;
        }

        currentNumClustersToAdd++;
        numClusterPairsTotal += currentCluster.getNeighbors().size();

        // If the thread is finished, write number of clusters and start new thread.
        if (numClusterPairsTotal >= numClusterPairsPerThread * (currentThread + 1)) {
          // Add the number of clusters for the finished thread.
          _clusterThreadPartition[currentThread].numClusters += currentNumClustersToAdd;
          currentNumClustersToAdd = 0;
          // Go to next thread!
          currentThread++;
          // if we are already at the end of all threads, go back to last thread!
          // this is a safety precaution and should not really matter.
          if (currentThread >= numThreads) {
            --currentThread;
            threadIsInitialized = true;
          } else {
            threadIsInitialized = false;
          }
        }
      }
    }
    if (not threadIsInitialized) {
      _clusterThreadPartition[currentThread] = {0, 0, 0};
    }
    // Make sure the last cluster range contains the rest of the clusters, even if there is not the perfect number left.
    if (currentNumClustersToAdd != 0) {
      _clusterThreadPartition[currentThread].numClusters += currentNumClustersToAdd;
    }
    // Theoretically, some threads may still remain. This ensures that their numClusters are set to 0.
    while (++currentThread < numThreads) {
      _clusterThreadPartition[currentThread] = {0, 0, 0};
    }
  }

  /**
   * If a particle is deleted, we want _isValid to be set to invalid, as the tower structure is invalidated.
   *
   * This function is not called, if a particle from the _particlesToAdd vector is deleted!
   */
  void notifyParticleDeleted() override {
    // this is potentially called from a threaded environment, so we have to make this atomic here!
    _isValid.store(ValidityState::invalid, std::memory_order::memory_order_relaxed);
  }

 private:
  /**
   * Helper function for the region iterators to determine bounds and towers to iterate over.
   * @param lowerCorner
   * @param higherCorner
   * @param behavior
   * @return
   */
  [[nodiscard]] auto getRegionIteratorHelper(const std::array<double, 3> &lowerCorner,
                                             const std::array<double, 3> &higherCorner,
                                             IteratorBehavior behavior) const {
    // Check all cells, as dummy particles are outside the domain they are only found if the search region is outside
    // the domain.
    const auto lowerCornerInBounds = utils::ArrayMath::max(lowerCorner, _haloBoxMin);
    const auto upperCornerInBounds = utils::ArrayMath::min(higherCorner, _haloBoxMax);

    if (not _builder) {
      // if no builder exists the clusters have not been built yet and all particles are stored in the first tower.
      return std::make_tuple(lowerCornerInBounds, upperCornerInBounds, std::vector<size_t>{0});
    }

    // Find towers intersecting the search region
    auto firstTowerCoords = _builder->getTowerCoordinates(lowerCornerInBounds);
    auto firstTowerIndex = _builder->towerIndex2DTo1D(firstTowerCoords[0], firstTowerCoords[1]);
    auto lastTowerCoords = _builder->getTowerCoordinates(upperCornerInBounds);

    std::array<size_t, 2> towersOfInterstPerDim;
    for (size_t dim = 0; dim < towersOfInterstPerDim.size(); ++dim) {
      // use ternary operators instead of abs because these are unsigned values
      towersOfInterstPerDim[dim] = firstTowerCoords[dim] > lastTowerCoords[dim]
                                       ? firstTowerCoords[dim] - lastTowerCoords[dim]
                                       : lastTowerCoords[dim] - firstTowerCoords[dim];
      // +1 because we want to include first AND last
      towersOfInterstPerDim[dim] += 1;
      // sanity check
      towersOfInterstPerDim[dim] = std::max(towersOfInterstPerDim[dim], static_cast<size_t>(1));
    }

    std::vector<size_t> towersOfInterest(towersOfInterstPerDim[0] * towersOfInterstPerDim[1]);

    auto towersOfInterestIterator = towersOfInterest.begin();
    for (size_t i = 0; i < towersOfInterstPerDim[1]; ++i) {
      std::iota(towersOfInterestIterator, towersOfInterestIterator + towersOfInterstPerDim[0],
                firstTowerIndex + i * _towersPerDim[0]);
      towersOfInterestIterator += towersOfInterstPerDim[0];
    }

    return std::make_tuple(lowerCornerInBounds, upperCornerInBounds, towersOfInterest);
  }

  /**
   * The number of particles in a full cluster.
   */
  size_t _clusterSize;

  /**
   * internal storage, particles are split into a grid in xy-dimension
   */
  std::vector<internal::ClusterTower<Particle>> _towers;

  /**
   * Dimensions of the 2D xy-grid.
   */
  std::array<size_t, 2> _towersPerDim{};

  /**
   * Side length of xy-grid.
   */
  double _towerSideLength{0.};
  double _towerSideLengthReciprocal{0.};

  /**
   * The number of clusters in the container.
   */
  size_t _numClusters;

  /**
   * The interaction length in number of towers it reaches.
   * static_cast<int>(std::ceil((this->getInteractionLength()) * _towerSideLengthReciprocal))
   */
  int _numTowersPerInteractionLength;

  /**
   * Contains all particles that should be added to the container during the next rebuild.
   */
  std::vector<Particle> _particlesToAdd;

  /**
   * Defines a partition of the clusters to a number of threads.
   */
  std::vector<ClusterRange> _clusterThreadPartition;

  /**
   * Minimum of the container.
   */
  std::array<double, 3> _boxMin{};

  /**
   * Maximum of the container.
   */
  std::array<double, 3> _boxMax{};

  /**
   * Minimum of the container including halo.
   */
  std::array<double, 3> _haloBoxMin{};

  /**
   * Maximum of the container including halo.
   */
  std::array<double, 3> _haloBoxMax{};

  /**
   * Cutoff.
   */
  double _cutoff{};

  /**
   * Skin.
   */
  double _skin{};

  /**
   * Enum to specify the validity of this container.
   */
  enum class ValidityState : unsigned char {
    invalid = 0,                 // nothing is valid.
    cellsValidListsInvalid = 1,  // only the cell structure is valid, but the lists are not.
    cellsAndListsValid = 2       // the cells and lists are valid
  };

  /**
   * Indicates, whether the current container structure (mainly for region iterators) and the verlet lists are valid.
   */
  std::atomic<ValidityState> _isValid{ValidityState::invalid};

  /**
   * The builder for the verlet cluster lists.
   */
  std::unique_ptr<internal::VerletClusterListsRebuilder<Particle>> _builder;

  /**
   * The flag manager of this container.
   */
  internal::UnknowingCellBorderAndFlagManager unknowingCellBorderAndFlagManager;
};

}  // namespace autopas
