/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>

#include "ClusterTowerBlock2D.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/NeighborListsBuffer.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/ParticleDeletedObserver.h"
#include "autopas/containers/cellTraversals/BalancedTraversal.h"
#include "autopas/containers/verletClusterLists/Cluster.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "autopas/containers/verletClusterLists/VerletClusterListsRebuilder.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas {

/**
 * Particles are divided into clusters.
 * The VerletClusterLists class uses neighborhood lists for each cluster to calculate pairwise interactions of
 * particles. It is optimized for a constant, i.e. particle independent, cutoff radius of the interaction.
 *
 * This Container does (currently?) not make use of the cellSizeFactor.
 *
 * The _particlesToAdd buffer structure is (currently) still necessary, even if the LogicHandler basically holds
 * the same buffer structure. In principle, moving the particles directly into one of the towers would be possible,
 * since particles are moved from LogicHandler to VCL only in a rebuild iteration. However, storing the particles in a
 * tower (e.g. tower0) is only possible very inefficiently, since several threads would write to this buffer at the same
 * time. Even if we could add the particles to tower0 efficiently, there are still problems in getRegionIterator(),
 * because this function is executed between the addition of particles and the actual rebuild. getRegionIterator()
 * expects that all particles are already sorted correctly into the towers (if we do not use _particlesToAdd).
 *
 * @note See VerletClusterListsRebuilder for the layout of the towers and clusters.
 *
 * @tparam Particle_T
 */
template <class Particle_T>
class VerletClusterLists : public ParticleContainerInterface<Particle_T>, public internal::ParticleDeletedObserver {
 public:
  /**
   * Type of the Particle.
   */
  using ParticleType = Particle_T;
  /**
   * Type of the ParticleCell which refers to the cluster towers.
   */
  using ParticleCellType = internal::ClusterTower<Particle_T>;

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
   * @param loadEstimator load estimation algorithm for balanced traversals.
   */
  VerletClusterLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff,
                     double skin, size_t clusterSize, LoadEstimatorOption loadEstimator = LoadEstimatorOption::none)
      : ParticleContainerInterface<Particle_T>(skin),
        _towerBlock{boxMin, boxMax, cutoff + skin},
        _clusterSize{clusterSize},
        _particlesToAdd(autopas_get_max_threads()),
        _cutoff{cutoff},
        _loadEstimator(loadEstimator) {
    // always have at least one tower.
    _towerBlock.addTower(_clusterSize);
  }

  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletClusterLists; }

  /**
   * Generates the load estimation function depending on _loadEstimator.
   * @return load estimator function object.
   */
  BalancedTraversal::EstimatorFunction getLoadEstimatorFunction() {
    // (Explicit) static cast required for Apple Clang (last tested version: 17.0.0)
    switch (static_cast<LoadEstimatorOption::Value>(this->_loadEstimator)) {
      case LoadEstimatorOption::neighborListLength: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          // the neighborListLength function defined for verletListsCells in not compatible with this container.
          unsigned long sum = 0;
          for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
            for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
              unsigned long cellLoad = 0;
              auto &tower = _towerBlock.getTowerByIndex2D(x, y);
              for (auto &cluster : tower.getClusters()) {
                if (cluster.getNeighbors()) {
                  cellLoad += cluster.getNeighbors()->size();
                }
              }
              sum += cellLoad;
            }
          }
          return sum;
        };
      }
      case LoadEstimatorOption::none:
        [[fallthrough]];
      default: {
        return
            [&](const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
                const std::array<unsigned long, 3> &upperCorner) { return 1; };
      }
    }
  }

  void computeInteractions(TraversalInterface *traversal) override {
    if (_isValid == ValidityState::cellsAndListsValid) {
      utils::ExceptionHandler::exception(
          "VerletClusterLists::computeInteractions(): Trying to do a pairwise iteration, even though verlet lists are "
          "not valid.");
    }
    auto *traversalInterface = dynamic_cast<VCLTraversalInterface<Particle_T> *>(traversal);
    if (traversalInterface) {
      traversalInterface->setClusterLists(*this);
      traversalInterface->setTowers(_towerBlock.getTowersRef());
    } else {
      utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in VerletClusterLists::computeInteractions. TraversalID: {}",
          traversal->getTraversalType());
    }
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }

    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
  }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    const auto particlesPerTower = (numParticles + numParticlesHaloEstimate) / _towerBlock.size();
    for (auto &tower : _towerBlock) {
      tower.reserve(particlesPerTower);
    }
  }

  /**
   * Adds the given particle to the container. rebuildTowersAndClusters() has to be called to have it actually sorted
   * in.
   * @param p The particle to add.
   */
  void addParticleImpl(const Particle_T &p) override {
    _isValid.store(ValidityState::invalid, std::memory_order_relaxed);
    _particlesToAdd[autopas_get_thread_num()].push_back(p);
  }

  void addHaloParticleImpl(const Particle_T &haloParticle) override {
    _isValid.store(ValidityState::invalid, std::memory_order_relaxed);
    _particlesToAdd[autopas_get_thread_num()].push_back(haloParticle);
  }

  bool updateHaloParticle(const Particle_T &haloParticle) override {
    using namespace autopas::utils::ArrayMath::literals;

    const auto &haloPos = haloParticle.getR();
    // this might be called from a parallel region so force this iterator to be sequential
    for (auto it = getRegionIterator(haloPos - (this->getVerletSkin() / 2.), haloPos + (this->getVerletSkin() / 2.),
                                     IteratorBehavior::halo | IteratorBehavior::forceSequential, nullptr);
         it.isValid(); ++it) {
      if (haloParticle.getID() == it->getID()) {
        // don't simply copy haloParticle over iter. This would trigger a dataRace with other regionIterators that
        // overlap with this region.
        it->setR(haloPos);
        it->setV(haloParticle.getV());
        it->setF(haloParticle.getF());
        return true;
      }
    }
    return false;
  }

  void deleteHaloParticles() override {
    // Step 1: Remove particles from _particlesToAdd
    for (auto &particleVec : _particlesToAdd) {
      for (size_t j = 0; j < particleVec.size();) {
        if (particleVec[j].isHalo()) {
          particleVec[j] = particleVec[particleVec.size() - 1];
          particleVec.pop_back();
        } else {
          ++j;
        }
      }
    }
    // Step 2: Remove particles from _towers
    bool deletedSomething = false;
    AUTOPAS_OPENMP(parallel for reduction(|| : deletedSomething) num_threads(autopas::autopas_get_preferred_num_threads()))
    // Thanks to clang 13 this has to be a index based loop instead of a range based
    for (size_t i = 0; i < _towerBlock.size(); ++i) {
      auto &tower = _towerBlock[i];
      const auto towerSize = tower.size();
      auto numTailDummies = tower.getNumTailDummyParticles();
      // iterate over all non-tail dummies. Avoid underflows.
      for (size_t j = 0; numTailDummies < towerSize and j < towerSize - numTailDummies;) {
        if (tower[j].isHalo()) {
          // swap-"delete"
          tower[j] = tower[towerSize - 1 - numTailDummies];
          // Since we can't pop the moved particle here mark it for deletion.
          internal::markParticleAsDeleted(tower[towerSize - 1 - numTailDummies]);
          ++numTailDummies;
          deletedSomething = true;
        } else {
          ++j;
        }
      }
      // if anything was marked for deletion actually delete it now.
      if (deletedSomething) {
        tower.deleteDummyParticles();
      }
    }
    if (deletedSomething) {
      _isValid.store(ValidityState::invalid, std::memory_order_relaxed);
    }
  }

  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             IteratorBehavior iteratorBehavior,
                                                             const std::array<double, 3> &boxMin,
                                                             const std::array<double, 3> &boxMax) const override {
    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }
  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             IteratorBehavior iteratorBehavior) const override {
    // this is not a region iter hence we stretch the bounding box to the numeric max
    constexpr std::array<double, 3> boxMin{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                                           std::numeric_limits<double>::lowest()};

    constexpr std::array<double, 3> boxMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                           std::numeric_limits<double>::max()};
    return getParticleImpl<false>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }

  /**
   * Container specific implementation for getParticle. See ParticleContainerInterface::getParticle().
   *
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @return tuple<ParticlePointer, CellIndex, ParticleIndex>
   */
  template <bool regionIter>
  std::tuple<const Particle_T *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
                                                                 IteratorBehavior iteratorBehavior,
                                                                 const std::array<double, 3> &boxMin,
                                                                 const std::array<double, 3> &boxMax) const {
    using namespace autopas::utils::ArrayMath::literals;

    // in this context cell == tower
    // catching the edge case that towers are not yet built -> Iterator jumps to additional vectors
    if (_towerBlock.empty()) {
      return {nullptr, 0, 0};
    }

    std::array<double, 3> boxMinWithSafetyMargin = boxMin;
    std::array<double, 3> boxMaxWithSafetyMargin = boxMax;
    if constexpr (regionIter) {
      // We extend the search box for cells here since particles might have moved
      boxMinWithSafetyMargin -= 0.5 * this->getVerletSkin();
      boxMaxWithSafetyMargin += 0.5 * this->getVerletSkin();
    }

    // first and last relevant cell index
    const auto [startCellIndex, endCellIndex] = [&]() -> std::tuple<size_t, size_t> {
      if constexpr (regionIter) {
        // We extend the search box for cells here since particles might have moved
        return {_towerBlock.getTowerIndex1DAtPosition(boxMinWithSafetyMargin),
                _towerBlock.getTowerIndex1DAtPosition(boxMaxWithSafetyMargin)};
      } else {
        if (not(iteratorBehavior & IteratorBehavior::halo)) {
          // only potentially owned region
          return {_towerBlock.getFirstOwnedTowerIndex(), _towerBlock.getLastOwnedTowerIndex()};
        } else {
          // whole range of cells
          return {0, _towerBlock.size() - 1};
        }
      }
    }();

    // if we are at the start of an iteration ...
    if (cellIndex == 0 and particleIndex == 0) {
      cellIndex =
          startCellIndex + ((iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num());
    }
    // abort if the start index is already out of bounds
    if (cellIndex >= _towerBlock.size()) {
      return {nullptr, 0, 0};
    }
    // check the data behind the indices
    if (particleIndex >= _towerBlock[cellIndex].getNumActualParticles() or
        not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            _towerBlock[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax)) {
      // either advance them to something interesting or invalidate them.
      std::tie(cellIndex, particleIndex) =
          advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax,
                                             boxMinWithSafetyMargin, boxMaxWithSafetyMargin, endCellIndex);
    }
    // shortcut if the given index doesn't exist
    if (cellIndex > endCellIndex) {
      return {nullptr, 0, 0};
    }
    const Particle_T *retPtr = &_towerBlock[cellIndex][particleIndex];

    return {retPtr, cellIndex, particleIndex};
  }

  bool deleteParticle(Particle_T &particle) override {
    // This function doesn't actually delete anything as it would mess up the clusters' references.
    internal::markParticleAsDeleted(particle);
    return false;
  }

  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {
    // This function doesn't actually delete anything as it would mess up the clusters' references.
    internal::markParticleAsDeleted(_towerBlock[cellIndex][particleIndex]);
    return false;
  }

  [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
    // First delete all halo particles.
    this->deleteHaloParticles();
    // Delete dummy particles.
    AUTOPAS_OPENMP(parallel for num_threads(autopas::autopas_get_preferred_num_threads()))
    for (size_t i = 0ul; i < _towerBlock.size(); ++i) {
      _towerBlock[i].deleteDummyParticles();
    }
    // Delete dummy particles also from the _particlesToAdd vector
    AUTOPAS_OPENMP(parallel for num_threads(autopas::autopas_get_preferred_num_threads()))
    for (size_t i = 0ul; i < _particlesToAdd.size(); ++i) {
      _particlesToAdd[i].erase(std::remove_if(_particlesToAdd[i].begin(), _particlesToAdd[i].end(),
                                              [](const auto &p) { return p.isDummy(); }),
                               _particlesToAdd[i].end());
    }

    // next find invalid particles
    std::vector<Particle_T> invalidParticles;

    // custom openmp reduction to concatenate all local vectors to one at the end of a parallel region
    AUTOPAS_OPENMP(declare reduction(
        vecMergeParticle : std::vector<Particle_T> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())))
    AUTOPAS_OPENMP(parallel reduction(vecMergeParticle : invalidParticles) num_threads(autopas::autopas_get_preferred_num_threads())) {
      for (auto iter = this->begin(IteratorBehavior::owned); iter.isValid(); ++iter) {
        if (not utils::inBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
          invalidParticles.push_back(*iter);
          internal::deleteParticle(iter);
        }
      }
    }
    _isValid.store(ValidityState::invalid, std::memory_order_relaxed);
    return invalidParticles;
  }

  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    using namespace autopas::utils::ArrayMath::literals;
    // Here, the towers might not yet be built, hence do not use members like _towerBlock.getTowersPerDim
    const auto boxSizeWithHalo = this->getHaloBoxMax() - this->getHaloBoxMin();
    const auto [towerSideLength, towersPerDim] = _towerBlock.estimateOptimalGridSideLength(
        this->getNumberOfParticles(IteratorBehavior::ownedOrHalo), _clusterSize);
    const std::array<double, 3> towerSize = {towerSideLength[0], towerSideLength[1],
                                             this->getHaloBoxMax()[2] - this->getHaloBoxMin()[2]};
    const std::array<unsigned long, 3> towerDimensions = {towersPerDim[0], towersPerDim[1], 1};
    return TraversalSelectorInfo(towerDimensions, this->getInteractionLength(), towerSize, _clusterSize);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    // Note: particlesToAddEmpty() can only be called if the container status is not invalid. If the status is set to
    // invalid, we do writing operations on _particlesToAdd and can not read from it without race conditions.
    if (_isValid != ValidityState::invalid) {
      // we call particlesToAddEmpty() as a sanity check to ensure there are actually no particles in _particlesToAdd if
      // the status is not invalid
      if (not particlesToAddEmpty(autopas_get_thread_num())) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::begin(): Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers + additionalVectors
      // from LogicHandler.
      return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle_T, true, false>::ParticleVecType additionalVectorsToPass;
      appendBuffersHelper(additionalVectors, additionalVectorsToPass);
      return ContainerIterator<Particle_T, true, false>(*this, behavior, &additionalVectorsToPass);
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version.
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors =
          nullptr) const override {
    // Note: particlesToAddEmpty() can only be called if the container status is not invalid. If the status is set to
    // invalid, we do writing operations on _particlesToAdd and can not read from from it without race conditions.
    if (_isValid != ValidityState::invalid) {
      // we call particlesToAddEmpty() as a sanity check to ensire there are actually no particles in _particlesToAdd if
      // the status is not invalid
      if (not particlesToAddEmpty(autopas_get_thread_num())) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::begin() const: Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers + additionalVectors
      // from LogicHandler.
      return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType additionalVectorsToPass;
      appendBuffersHelper(additionalVectors, additionalVectorsToPass);
      return ContainerIterator<Particle_T, false, false>(*this, behavior, &additionalVectorsToPass);
    }
  }

  /**
   * @copydoc autopas::LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) {
    for (auto &tower : _towerBlock) {
      tower.forEach(forEachLambda, behavior);
    }
    for (auto &vector : this->_particlesToAdd) {
      for (auto &particle : vector) {
        if (behavior.contains(particle)) {
          forEachLambda(particle);
        }
      }
    }
  }

  /**
   * @copydoc autopas::LinkedCells::forEach()
   * @note const version.
   * @note This function additionally iterates over the _particlesToAdd vector if the tower-structure isn't valid.
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::forEach() const: Error: particle container is valid, but _particlesToAdd isn't "
            "empty!");
      }
    }

    // If the particles are sorted into the towers, we can simply use the iteration over towers.
    for (auto &tower : _towerBlock) {
      tower.forEach(forEachLambda, behavior);
    }

    // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
    if (_isValid == ValidityState::invalid) {
      for (auto &particlesToAddPerThread : _particlesToAdd) {
        for (auto &particle : particlesToAddPerThread) {
          if (behavior.contains(particle)) {
            forEachLambda(particle);
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::LinkedCells::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) {
    for (auto &tower : _towerBlock) {
      tower.reduce(reduceLambda, result, behavior);
    }
    for (auto &vector : this->_particlesToAdd) {
      for (auto &p : vector) {
        if (behavior.contains(p)) {
          reduceLambda(p, result);
        }
      }
    }
  }

  /**
   * @copydoc autopas::LinkedCells::reduce()
   * @note const version.
   * @note This function additionally iterates over the _particlesToAdd vector if the tower-structure isn't valid.
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result,
              IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::reduce() const: Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
    }

    // If the particles are sorted into the towers, we can simply use the iteration over towers.
    for (auto tower : _towerBlock) {
      tower.reduce(reduceLambda, result, behavior);
    }

    if (_isValid == ValidityState::invalid) {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      for (auto &particlesToAddPerThread : _particlesToAdd) {
        for (auto &particle : particlesToAddPerThread) {
          if (behavior.contains(particle)) {
            reduceLambda(particle, result);
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   */
  [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors = nullptr) override {
    // Note: particlesToAddEmpty() can only be called if the container status is not invalid. If the status is set to
    // invalid, we do writing operations on _particlesToAdd and can not read from from it without race conditions.
    if (_isValid != ValidityState::invalid) {
      // we call particlesToAddEmpty() as a sanity check to ensure there are actually no particles in _particlesToAdd
      // if the status is not invalid
      if (not particlesToAddEmpty(autopas_get_thread_num())) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::reduce() const: Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers + additionalVectors
      // from LogicHandler.
      return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle_T, true, true>::ParticleVecType additionalVectorsToPass;
      appendBuffersHelper(additionalVectors, additionalVectorsToPass);
      return ContainerIterator<Particle_T, true, true>(*this, behavior, &additionalVectorsToPass, lowerCorner,
                                                       higherCorner);
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version.
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors =
          nullptr) const override {
    // Note: particlesToAddEmpty() can only be called if the container status is not invalid. If the status is set to
    // invalid, we do writing operations on _particlesToAdd and can not read from from it without race conditions.
    if (_isValid != ValidityState::invalid) {
      // we call particlesToAddEmpty() as a sanity check to ensire there are actually no particles in _particlesToAdd if
      // the status is not invalid
      if (not particlesToAddEmpty(autopas_get_thread_num())) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::getRegionIterator() const: Error: particle container is valid, but _particlesToAdd "
            "isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers + additionalVectors
      // from LogicHandler.
      return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle_T, false, true>::ParticleVecType additionalVectorsToPass;
      appendBuffersHelper(additionalVectors, additionalVectorsToPass);
      return ContainerIterator<Particle_T, false, true>(*this, behavior, &additionalVectorsToPass, lowerCorner,
                                                        higherCorner);
    }
  }

  /**
   * @copydoc autopas::LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner,
                       IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) {
    for (size_t i = 0; i < _towerBlock.size(); ++i) {
      if (_towerBlock.ignoreCellForIteration(i, behavior)) {
        continue;
      }
      auto &tower = _towerBlock[i];
      const auto [towerLowCorner, towerHighCorner] = _towerBlock.getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      }
    }
    for (auto &vector : _particlesToAdd) {
      for (auto &particle : vector) {
        if (behavior.contains(particle)) {
          if (utils::inBox(particle.getR(), lowerCorner, higherCorner)) {
            forEachLambda(particle);
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::LinkedCells::forEachInRegion()
   * @note const version.
   * @note This function additionally iterates over the _particlesToAdd vector if the tower-structure isn't valid.
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner,
                       IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::forEachInRegion() const: Error: particle container is valid, but _particlesToAdd "
            "isn't empty!");
      }
    }

    // If the particles are sorted into the towers, we can simply use the iteration over towers.
    for (size_t i = 0; i < _towerBlock.size(); ++i) {
      if (_towerBlock.ignoreCellForIteration(i, behavior)) {
        continue;
      }
      auto &tower = _towerBlock[i];
      const auto [towerLowCorner, towerHighCorner] = _towerBlock.getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      }
    }

    if (_isValid == ValidityState::invalid) {
      // If the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      for (auto &particlesToAddPerThread : _particlesToAdd) {
        for (auto &particle : particlesToAddPerThread) {
          if (behavior.contains(particle)) {
            if (utils::inBox(particle.getR(), lowerCorner, higherCorner)) {
              forEachLambda(particle);
            }
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner,
                      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) {
    for (size_t i = 0; i < _towerBlock.size(); ++i) {
      if (_towerBlock.ignoreCellForIteration(i, behavior)) {
        continue;
      }
      auto &tower = _towerBlock[i];
      const auto [towerLowCorner, towerHighCorner] = _towerBlock.getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
      }
    }
    for (auto &vector : _particlesToAdd) {
      for (auto &particle : vector) {
        if (behavior.contains(particle)) {
          if (utils::inBox(particle.getR(), lowerCorner, higherCorner)) {
            reduceLambda(particle, result);
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::LinkedCells::reduceInRegion()
   * @note const version.
   * @note This function additionally iterates over the _particlesToAdd vector if the tower-structure isn't valid.
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner,
                      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::reduceInRegion() const: Error: particle container is valid, but _particlesToAdd isn't "
            "empty!");
      }
    }
    // If the particles are sorted into the towers, we can simply use the iteration over towers.
    for (size_t i = 0; i < _towerBlock.size(); ++i) {
      if (_towerBlock.ignoreCellForIteration(i, behavior)) {
        continue;
      }
      auto &tower = _towerBlock[i];
      const auto [towerLowCorner, towerHighCorner] = _towerBlock.getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
      }
    }

    if (_isValid == ValidityState::invalid) {
      // If the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      for (auto &particlesToAddPerThread : _particlesToAdd) {
        for (auto &particle : particlesToAddPerThread) {
          if (behavior.contains(particle)) {
            if (utils::inBox(particle.getR(), lowerCorner, higherCorner)) {
              reduceLambda(particle, result);
            }
          }
        }
      }
    }
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // The builder might have a different newton3 choice than the traversal. This typically only happens in unit tests
    // when rebuildTowersAndClusters() was not called explicitly.
    if (_isValid == ValidityState::invalid or traversal->getUseNewton3() != _builder->getNewton3()) {
      // clear the lists buffer because clusters will be recreated
      _neighborLists.clear();
      rebuildTowersAndClusters(traversal->getUseNewton3());
    }
    _builder->rebuildNeighborListsAndFillClusters();

    auto *clusterTraversalInterface = dynamic_cast<VCLTraversalInterface<Particle_T> *>(traversal);
    if (clusterTraversalInterface) {
      if (clusterTraversalInterface->needsStaticClusterThreadPartition()) {
        calculateClusterThreadPartition();
      }
    } else {
      utils::ExceptionHandler::exception(
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
   * Get the number of all particles stored in this container (owned + halo + dummy).
   * @return number of particles stored in this container (owned + halo + dummy).
   */
  [[nodiscard]] size_t size() const override {
    size_t sum = std::accumulate(_towerBlock.begin(), _towerBlock.end(), 0,
                                 [](size_t acc, const auto &tower) { return acc + tower.size(); });
    sum = std::accumulate(_particlesToAdd.begin(), _particlesToAdd.end(), sum,
                          [](size_t acc, const auto &buffer) { return acc + buffer.size(); });
    return sum;
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getNumberOfParticles()
   */
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior) const override {
    // sum up all particles in towers that fulfill behavior
    size_t sum = std::accumulate(_towerBlock.begin(), _towerBlock.end(), 0, [&behavior](size_t acc, const auto &tower) {
      return acc + tower.getNumberOfParticles(behavior);
    });

    // Since we can not directly insert particles into towers without a rebuild of the whole data structure,
    // _particlesToAdd is used to store all these particles temporarily until the next rebuild inserts them into the
    // towers data structure. However, these particles already belong to the respective tower, so we have to count them
    // as well.
    sum = std::accumulate(
        _particlesToAdd.begin(), _particlesToAdd.end(), sum, [&behavior](size_t acc, const auto &buffer) {
          return acc +
                 (std::count_if(buffer.begin(), buffer.end(), [&behavior](auto p) { return behavior.contains(p); }));
        });

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
  auto getTowerSideLength() const { return _towerBlock.getTowerSideLength(); }

  /**
   * Returns the number of grids per dimension on the container.
   * @return the number of grids per dimension on the container.
   */
  auto getTowersPerDimension() const { return _towerBlock.getTowersPerDim(); }

  /**
   * Returns the number of particles in each cluster.
   * @return the number of particles in each cluster.
   */
  auto getClusterSize() const { return _clusterSize; }

  /**
   * Returns the towers per interaction length. That is how many towers fit into one interaction length rounded up.
   * @return the number of towers per interaction length.
   */
  auto getNumTowersPerInteractionLength() const { return _towerBlock.getNumTowersPerInteractionLength(); }

  /**
   * Loads all particles of the container in their correct SoA and generates the SoAViews for the clusters.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for loading the particles into the SoA.
   */
  template <class Functor>
  void loadParticlesIntoSoAs(Functor *functor) {
    const auto numTowers = _towerBlock.size();
    /// @todo: find sensible chunksize
    AUTOPAS_OPENMP(parallel for schedule(dynamic) num_threads(autopas::autopas_get_preferred_num_threads()))
    for (size_t index = 0; index < numTowers; index++) {
      _towerBlock[index].loadSoA(functor);
    }
  }

  /**
   * Extracts all SoAs of the container into the particles.
   * @tparam Functor The type of the functor to use.
   * @param functor The functor to use for extracting the SoAs into the particles..
   */
  template <class Functor>
  void extractParticlesFromSoAs(Functor *functor) {
    const auto numTowers = _towerBlock.size();
    /// @todo: find sensible chunksize
    AUTOPAS_OPENMP(parallel for schedule(dynamic) num_threads(autopas::autopas_get_preferred_num_threads()))
    for (size_t index = 0; index < numTowers; index++) {
      _towerBlock[index].extractSoA(functor);
    }
  }

  /**
   * Returns a reference to the tower for the given tower grid coordinates.
   * @param x The x-th tower in x direction.
   * @param y The y-th tower in y direction.
   * @return a reference to the tower for the given tower grid coordinates.
   */
  internal::ClusterTower<Particle_T> &getTowerByIndex(size_t x, size_t y) {
    return _towerBlock.getTowerByIndex2D(x, y);
  }

  /**
   * Getter for the cell block.
   *
   * @note This is only used for testing.
   *
   * @return
   */
  internal::ClusterTowerBlock2D<Particle_T> &getTowerBlock() { return _towerBlock; }

  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override { return _towerBlock.getBoxMax(); }

  /**
   * Get the upper corner of the halo box.
   * @return the upper corner of the halo box.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMax() const { return _towerBlock.getHaloBoxMax(); }

  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override { return _towerBlock.getBoxMin(); }

  /**
   * Get the lower corner of the halo box.
   * @return the lower corner of the halo box.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMin() const { return _towerBlock.getHaloBoxMin(); }

  [[nodiscard]] double getCutoff() const override { return _cutoff; }

  void setCutoff(double cutoff) override { _cutoff = cutoff; }

  [[nodiscard]] double getVerletSkin() const override { return this->_skin; }

  /**
   * Set the verlet skin length for the container.
   * @param skin
   */
  void setSkin(double skin) { this->_skin = skin; }

  [[nodiscard]] double getInteractionLength() const override { return _cutoff + this->getVerletSkin(); }

  void deleteAllParticles() override {
    _isValid.store(ValidityState::invalid, std::memory_order_relaxed);
    std::for_each(_particlesToAdd.begin(), _particlesToAdd.end(), [](auto &buffer) { buffer.clear(); });
    std::for_each(_towerBlock.begin(), _towerBlock.end(), [](auto &tower) { tower.clear(); });
  }

  /**
   * Get the neighbor lists buffer object.
   * @return
   */
  const typename internal::VerletClusterListsRebuilder<Particle_T>::NeighborListsBuffer_T &getNeighborLists() const {
    return _neighborLists;
  }

  /**
   * Initializes a new VerletClusterListsRebuilder and uses it to rebuild the towers and the clusters.
   * This function sets the container structure to valid.
   * @param newton3 Indicate whether the VerletClusterRebuilder should consider newton3 or not.
   */
  void rebuildTowersAndClusters(bool newton3) {
    using namespace utils::ArrayMath::literals;
    // collect all particles to add from across the thread buffers
    typename decltype(_particlesToAdd)::value_type particlesToAdd;
    const size_t numParticlesToAdd =
        std::accumulate(_particlesToAdd.begin(), _particlesToAdd.end(), 0,
                        [](size_t acc, const auto &buffer) { return acc + buffer.size(); });
    particlesToAdd.reserve(numParticlesToAdd);
    std::for_each(_particlesToAdd.begin(), _particlesToAdd.end(), [&](auto &particlesBuffer) {
      particlesToAdd.insert(particlesToAdd.end(), particlesBuffer.begin(), particlesBuffer.end());
      particlesBuffer.clear();
    });

    const double interactionLength = _cutoff + this->_skin;
    _builder = std::make_unique<internal::VerletClusterListsRebuilder<Particle_T>>(
        _towerBlock, particlesToAdd, _neighborLists, _clusterSize, interactionLength * interactionLength, newton3);

    _numClusters = _builder->rebuildTowersAndClusters();

    _isValid.store(ValidityState::cellsValidListsInvalid, std::memory_order_relaxed);
    for (auto &tower : _towerBlock) {
      tower.setParticleDeletionObserver(this);
    }
  }

  /**
   * Helper method to sequentially iterate over all owned clusters.
   * @tparam LoopBody The type of the lambda to execute for all clusters.
   * @param loopBody The lambda to execute for all clusters. Parameters given is internal::Cluster& cluster.
   */
  template <class LoopBody>
  void traverseClustersSequential(LoopBody &&loopBody) {
    for (size_t x = 0; x < _towerBlock.getTowersPerDim()[0]; x++) {
      for (size_t y = 0; y < _towerBlock.getTowersPerDim()[1]; y++) {
        auto &tower = _towerBlock.getTowerByIndex2D(x, y);
        for (auto clusterIter = tower.getFirstOwnedCluster(); clusterIter < tower.getFirstTailHaloCluster();
             ++clusterIter) {
          loopBody(*clusterIter);
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
    const auto towersPerDimX = _towerBlock.getTowersPerDim()[0];
    const auto towersPerDimY = _towerBlock.getTowersPerDim()[1];
    /// @todo: find sensible chunksize
    AUTOPAS_OPENMP(parallel for schedule(dynamic) collapse(2) num_threads(autopas::autopas_get_preferred_num_threads()))
    for (size_t x = 0; x < towersPerDimX; x++) {
      for (size_t y = 0; y < towersPerDimY; y++) {
        auto &tower = _towerBlock.getTowerByIndex2D(x, y);

        for (auto clusterIter = tower.getFirstOwnedCluster(); clusterIter < tower.getFirstTailHaloCluster();
             ++clusterIter) {
          loopBody(*clusterIter);
        }
      }
    }
  }

 protected:
  /**
   * Calculates a cluster thread partition that aims to give each thread about the same amount of cluster pair
   * interactions, if each thread handles the neighbors of all clusters it gets assigned.
   */
  void calculateClusterThreadPartition() {
    size_t numClusterPairs = 0;
    this->template traverseClusters<false>(
        [&numClusterPairs](auto &cluster) { numClusterPairs += cluster.getNeighbors()->size(); });

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
    for (size_t currentTowerIndex = 0; currentTowerIndex < _towerBlock.size(); currentTowerIndex++) {
      auto &currentTower = _towerBlock[currentTowerIndex];
      const auto firstOwnedClusterIndex =
          std::distance(currentTower.getClusters().begin(), currentTower.getFirstOwnedCluster());
      const auto firstTailHaloClusterIndex =
          std::distance(currentTower.getClusters().begin(), currentTower.getFirstTailHaloCluster());
      for (size_t currentClusterInTower = firstOwnedClusterIndex; currentClusterInTower < firstTailHaloClusterIndex;
           ++currentClusterInTower) {
        auto &currentCluster = currentTower.getCluster(currentClusterInTower);

        // If on a new thread, start with the clusters for this thread here.
        if (not threadIsInitialized) {
          _clusterThreadPartition[currentThread] = {currentTowerIndex, currentClusterInTower, 0};
          threadIsInitialized = true;
        }

        currentNumClustersToAdd++;
        numClusterPairsTotal += currentCluster.getNeighbors()->size();

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
    _isValid.store(ValidityState::invalid, std::memory_order_relaxed);
  }

 private:
  /**
   * Given a pair of cell-/particleIndex and iterator restrictions either returns the next indices that match these
   * restrictions or indices that are out of bounds (e.g. cellIndex >= cells.size())
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin The actual search box min
   * @param boxMax The actual search box max
   * @param boxMinWithSafetyMargin Search box min that includes a surrounding of skin
   * @param boxMaxWithSafetyMargin Search box max that includes a surrounding of skin
   * @param endCellIndex
   * @return tuple<cellIndex, particleIndex>
   */
  template <bool regionIter>
  [[nodiscard]] std::tuple<size_t, size_t> advanceIteratorIndices(
      size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin,
      const std::array<double, 3> &boxMax, const std::array<double, 3> &boxMinWithSafetyMargin,
      const std::array<double, 3> &boxMaxWithSafetyMargin, size_t endCellIndex) const {
    // Finding the indices for the next particle
    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_preferred_num_threads();

    // helper function to determine if the cell can even contain particles of interest to the iterator
    auto towerIsRelevant = [&]() -> bool {
      // special case: Towers are not yet built, then the tower acting as buffer is always relevant.
      if (_towerBlock.size() == 1) {
        return true;
      }
      bool isRelevant = true;
      if constexpr (regionIter) {
        // is the cell in the region?
        const auto [towerLowCorner, towerHighCorner] = _towerBlock.getTowerBoundingBox(cellIndex);
        isRelevant =
            utils::boxesOverlap(towerLowCorner, towerHighCorner, boxMinWithSafetyMargin, boxMaxWithSafetyMargin);
      }
      return isRelevant;
    };

    do {
      // advance to the next particle
      ++particleIndex;
      // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.

      // If cell has wrong type, or there are no more particles in this cell jump to the next
      while (not towerIsRelevant() or particleIndex >= _towerBlock[cellIndex].getNumActualParticles()) {
        cellIndex += stride;
        particleIndex = 0;

        // If we notice that there is nothing else to look at set invalid values, so we get a nullptr next time and
        // break.
        if (cellIndex > endCellIndex) {
          return {std::numeric_limits<size_t>::max(), particleIndex};
        }
      }
    } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
        _towerBlock[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax));

    // the indices returned at this point should always be valid
    return {cellIndex, particleIndex};
  }

  internal::ClusterTowerBlock2D<Particle_T> _towerBlock;

  /**
   * load estimation algorithm for balanced traversals.
   */
  autopas::LoadEstimatorOption _loadEstimator;

  /**
   * The number of particles in a full cluster.
   */
  size_t _clusterSize;

  /**
   * The number of clusters in the container.
   */
  size_t _numClusters{0};

  /**
   * Contains all particles that should be added to the container during the next rebuild.
   * Outer vector is for Thread buffer to allow parallel particle insertion.
   * This has to be a mutable so we can call appendBuffersHelper() from const and non-const functions.
   */
  mutable std::vector<std::vector<Particle_T>> _particlesToAdd;

  /**
   * Checks if there are particles in the buffers of _particlesToAdd.
   * @param bufferID the buffer ID to check for emptiness. If bufferID == -1, all buffers are checked
   * @return true if all buffers are empty or if one of the specified buffers is empty.
   */
  [[nodiscard]] bool particlesToAddEmpty(int bufferID = -1) const {
    if (bufferID == -1) {
      for (auto &threadBuffer : _particlesToAdd) {
        if (not threadBuffer.empty()) {
          return false;
        }
      }
      return true;
    } else {
      return _particlesToAdd[bufferID].empty();
    }
  }

  /**
   * Helper function for begin() and getRegionIterator() that merges all buffers from _particlesToAdd and
   * additionalVectors into a single buffer
   *
   * @tparam VecVec Type of datastructure of additional vectors. Expected is a vector of vectors.
   * @param additionalVectors Additional vectors from LogicHandler.
   * @param outVec The buffer where additionalVectors + _particlesToAdd will be stored.
   */
  template <class VecVec>
  void appendBuffersHelper(VecVec *additionalVectors, VecVec &outVec) const {
    if (additionalVectors) {
      outVec.reserve(_particlesToAdd.size() + additionalVectors->size());
      outVec.insert(outVec.end(), additionalVectors->begin(), additionalVectors->end());
    } else {
      outVec.reserve(_particlesToAdd.size());
    }
    for (auto &vec : _particlesToAdd) {
      outVec.push_back(&vec);
    }
  }

  /**
   * Defines a partition of the clusters to a number of threads.
   */
  std::vector<ClusterRange> _clusterThreadPartition;

  /**
   * Cutoff.
   */
  double _cutoff{};

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
  std::unique_ptr<internal::VerletClusterListsRebuilder<Particle_T>> _builder;

  /**
   * Structure to provide persistent memory for neighbor lists. Will be filled by the builder.
   */
  typename internal::VerletClusterListsRebuilder<Particle_T>::NeighborListsBuffer_T _neighborLists{};
};

}  // namespace autopas
