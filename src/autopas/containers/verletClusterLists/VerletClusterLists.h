/**
 * @file VerletClusterLists.h
 * @author nguyen
 * @date 14.10.18
 */

#pragma once

#include <algorithm>
#include <cmath>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/ParticleDeletedObserver.h"
#include "autopas/containers/cellPairTraversals/BalancedTraversal.h"
#include "autopas/containers/verletClusterLists/ClusterTower.h"
#include "autopas/containers/verletClusterLists/VerletClusterListsRebuilder.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Timer.h"
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
 * In the current solution, we can make use of the fact that _particlesToAdd only contains particles in a
 * rebuild-iteration and the additionalVectors only contain particles in a non-rebuild-iteration. This means that we can
 * save expensive buffer allocations in begin() and getRegionIterator() in valid container-status iterations, by just
 * passing the additionalVectors instead of concatenating them with _particlesToAdd.
 *
 * @note See VerletClusterListsRebuilder for the layout of the towers and clusters.
 *
 * @tparam Particle
 */
template <class Particle>
class VerletClusterLists : public ParticleContainerInterface<Particle>, public internal::ParticleDeletedObserver {
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
   * @param skinPerTimestep The skin radius per Timestep.
   * @param rebuildFrequency The rebuild Frequency.
   * @param clusterSize Number of particles per cluster.
   * @param loadEstimator load estimation algorithm for balanced traversals.
   */
  VerletClusterLists(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double cutoff,
                     double skinPerTimestep, unsigned int rebuildFrequency, size_t clusterSize,
                     LoadEstimatorOption loadEstimator = LoadEstimatorOption::none)
      : ParticleContainerInterface<Particle>(),
        _clusterSize{clusterSize},
        _numClusters{0},
        _numTowersPerInteractionLength{0},
        _particlesToAdd(autopas_get_max_threads()),
        _boxMin{boxMin},
        _boxMax{boxMax},
        _haloBoxMin{utils::ArrayMath::subScalar(boxMin, cutoff + skinPerTimestep * rebuildFrequency)},
        _haloBoxMax{utils::ArrayMath::addScalar(boxMax, cutoff + skinPerTimestep * rebuildFrequency)},
        _cutoff{cutoff},
        _skinPerTimestep{skinPerTimestep},
        _rebuildFrequency{rebuildFrequency},
        _loadEstimator(loadEstimator) {
    // always have at least one tower.
    _towers.push_back(internal::ClusterTower<Particle>(_clusterSize));
  }

  CellType getParticleCellTypeEnum() override { return CellType::ClusterTower; };

  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::verletClusterLists; }

  /**
   * Generates the load estimation function depending on _loadEstimator.
   * @return load estimator function object.
   */
  BalancedTraversal::EstimatorFunction getLoadEstimatorFunction() {
    switch (this->_loadEstimator) {
      case LoadEstimatorOption::neighborListLength: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          // the neighborListLength function defined for verletListsCells in not compatible with this container.
          unsigned long sum = 0;
          for (unsigned long x = lowerCorner[0]; x <= upperCorner[0]; x++) {
            for (unsigned long y = lowerCorner[1]; y <= upperCorner[1]; y++) {
              unsigned long cellLoad = 0;
              auto &tower = getTowerByIndex(x, y);
              for (auto &cluster : tower.getClusters()) {
                cellLoad += cluster.getNeighbors().size();
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
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    const auto particlesPerTower = (numParticles + numParticlesHaloEstimate) / _towers.size();
    for (auto &tower : _towers) {
      tower.reserve(particlesPerTower);
    }
  }

  /**
   * Adds the given particle to the container. rebuildTowersAndClusters() has to be called to have it actually sorted
   * in.
   * @param p The particle to add.
   */
  void addParticleImpl(const Particle &p) override {
    _isValid.store(ValidityState::invalid, std::memory_order::memory_order_relaxed);
    _particlesToAdd[autopas_get_thread_num()].push_back(p);
  }

  void addHaloParticleImpl(const Particle &haloParticle) override {
    _isValid.store(ValidityState::invalid, std::memory_order::memory_order_relaxed);
    Particle copy = haloParticle;
    copy.setOwnershipState(OwnershipState::halo);
    _particlesToAdd[autopas_get_thread_num()].push_back(copy);
  }

  bool updateHaloParticle(const Particle &haloParticle) override {
    using namespace autopas::utils::ArrayMath::literals;

    Particle pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);

    typename ContainerIterator<Particle, true, true>::ParticleVecType additionalVectors;
    additionalVectors.reserve(_particlesToAdd.size());
    for (auto &v : _particlesToAdd) {
      additionalVectors.push_back(&v);
    }

    // this might be called from a parallel region so force this iterator to be sequential
    for (auto it =
             getRegionIterator(pCopy.getR() - (this->getVerletSkin() / 2), pCopy.getR() + (this->getVerletSkin() / 2),
                               IteratorBehavior::halo | IteratorBehavior::forceSequential, &additionalVectors);
         it.isValid(); ++it) {
      if (pCopy.getID() == it->getID()) {
        *it = pCopy;
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
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for reduction(|| : deletedSomething)
#endif
    for (size_t i = 0; i < _towers.size(); ++i) {
      auto &tower = _towers[i];
      const auto towerSize = tower.getNumAllParticles();
      auto numTailDummies = tower.getNumTailDummyParticles();
      // iterate over all non-tail dummies.
      for (size_t j = 0; j < towerSize - numTailDummies;) {
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
      _isValid.store(ValidityState::invalid, std::memory_order::memory_order_relaxed);
    }
  }

  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                           IteratorBehavior iteratorBehavior,
                                                           const std::array<double, 3> &boxMin,
                                                           const std::array<double, 3> &boxMax) const override {
    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
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
  std::tuple<const Particle *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
                                                               IteratorBehavior iteratorBehavior,
                                                               const std::array<double, 3> &boxMin,
                                                               const std::array<double, 3> &boxMax) const {
    // in this context cell == tower
    // catching the edge case that towers are not yet built -> Iterator jumps to additional vectors
    if (_towers.empty()) {
      return {nullptr, 0, 0};
    }

    // first and last relevant cell index
    const auto [startCellIndex, endCellIndex] = [&]() -> std::tuple<size_t, size_t> {
      if constexpr (regionIter) {
        return {getTowerIndex1DAtPosition(boxMin), getTowerIndex1DAtPosition(boxMax)};
      } else {
        // whole range of cells
        return {0, this->_towers.size() - 1};
      }
    }();

    // if we are at the start of an iteration ...
    if (cellIndex == 0 and particleIndex == 0) {
      cellIndex =
          startCellIndex + ((iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num());
    }
    // abort if the start index is already out of bounds
    if (cellIndex >= this->_towers.size()) {
      return {nullptr, 0, 0};
    }
    // check the data behind the indices
    if (particleIndex >= this->_towers[cellIndex].numParticles() or
        not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            this->_towers[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax)) {
      // either advance them to something interesting or invalidate them.
      std::tie(cellIndex, particleIndex) =
          advanceIteratorIndices<regionIter>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax, endCellIndex);
    }

    // shortcut if the given index doesn't exist
    if (cellIndex > endCellIndex) {
      return {nullptr, 0, 0};
    }
    const Particle *retPtr = &this->_towers[cellIndex][particleIndex];

    return {retPtr, cellIndex, particleIndex};
  }

  bool deleteParticle(Particle &particle) override {
    // This function doesn't actually delete anything as it would mess up the clusters' references.
    internal::markParticleAsDeleted(particle);
    return false;
  }

  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {
    // This function doesn't actually delete anything as it would mess up the clusters' references.
    internal::markParticleAsDeleted(_towers[cellIndex][particleIndex]);
    return false;
  }

  [[nodiscard]] std::vector<Particle> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
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
#pragma omp declare reduction(vecMergeParticle : std::vector<Particle> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel reduction(vecMergeParticle : invalidParticles)
#endif
    {
      for (auto iter = this->begin(IteratorBehavior::owned); iter.isValid(); ++iter) {
        if (not utils::inBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
          invalidParticles.push_back(*iter);
          internal::deleteParticle(iter);
        }
      }
    }
    _isValid.store(ValidityState::invalid, std::memory_order::memory_order_relaxed);
    return invalidParticles;
  }

  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    using namespace autopas::utils::ArrayMath::literals;

    auto boxSizeWithHalo = this->getHaloBoxMax() - this->getHaloBoxMin();
    auto towerSideLength = internal::VerletClusterListsRebuilder<Particle>::estimateOptimalGridSideLength(
        this->getNumberOfParticles(), boxSizeWithHalo, _clusterSize);
    auto towersPerDim =
        internal::VerletClusterListsRebuilder<Particle>::calculateTowersPerDim(boxSizeWithHalo, 1.0 / towerSideLength);
    const std::array<double, 3> towerSize = {towerSideLength, towerSideLength,

                                             this->getHaloBoxMax()[2] - this->getHaloBoxMin()[2]};
    const std::array<unsigned long, 3> towerDimensions = {towersPerDim[0], towersPerDim[1], 1};
    return TraversalSelectorInfo(towerDimensions, this->getInteractionLength(), towerSize, _clusterSize);
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   */
  [[nodiscard]] ContainerIterator<Particle, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::begin(): Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers.
      return ContainerIterator<Particle, true, false>(*this, behavior, additionalVectors);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle, true, false>::ParticleVecType additionalVectorsToPass;
      if (not additionalVectorsEmpty<true, false>(additionalVectors)) {
        additionalVectorsToPass.reserve(_particlesToAdd.size() + additionalVectors->size());
        additionalVectorsToPass.insert(additionalVectorsToPass.end(), additionalVectors->begin(),
                                       additionalVectors->end());
      } else {
        additionalVectorsToPass.reserve(_particlesToAdd.size());
      }
      for (auto &vec : _particlesToAdd) {
        additionalVectorsToPass.push_back(&vec);
      }
      return ContainerIterator<Particle, true, false>(*this, behavior, &additionalVectorsToPass);
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version.
   */
  [[nodiscard]] ContainerIterator<Particle, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle, false, false>::ParticleVecType *additionalVectors = nullptr) const override {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::begin() const: Error: particle container is valid, but _particlesToAdd isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers.
      return ContainerIterator<Particle, false, false>(*this, behavior, additionalVectors);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle, false, false>::ParticleVecType additionalVectorsToPass;
      if (not additionalVectorsEmpty<false, false>(additionalVectors)) {
        additionalVectorsToPass.reserve(_particlesToAdd.size() + additionalVectors->size());
        additionalVectorsToPass.insert(additionalVectorsToPass.end(), additionalVectors->begin(),
                                       additionalVectors->end());
      } else {
        additionalVectorsToPass.reserve(_particlesToAdd.size());
      }
      for (auto &vec : _particlesToAdd) {
        additionalVectorsToPass.push_back(&vec);
      }
      return ContainerIterator<Particle, false, false>(*this, behavior, &additionalVectorsToPass);
    }
  }

  /**
   * @copydoc autopas::LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) {
    for (auto &tower : this->_towers) {
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
    for (auto &tower : this->_towers) {
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
    for (auto &tower : this->_towers) {
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
    for (auto tower : this->_towers) {
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
  [[nodiscard]] ContainerIterator<Particle, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle, true, true>::ParticleVecType *additionalVectors) override {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::getRegionIterator(): Error: particle container is valid, but _particlesToAdd "
            "isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers.
      return ContainerIterator<Particle, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle, true, true>::ParticleVecType additionalVectorsToPass;
      if (not additionalVectorsEmpty<true, true>(additionalVectors)) {
        additionalVectorsToPass.reserve(_particlesToAdd.size() + additionalVectors->size());
        additionalVectorsToPass.insert(additionalVectorsToPass.end(), additionalVectors->begin(),
                                       additionalVectors->end());
      } else {
        additionalVectorsToPass.reserve(_particlesToAdd.size());
      }
      for (auto &vec : _particlesToAdd) {
        additionalVectorsToPass.push_back(&vec);
      }
      return ContainerIterator<Particle, true, true>(*this, behavior, &additionalVectorsToPass, lowerCorner,
                                                     higherCorner);
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version.
   */
  [[nodiscard]] ContainerIterator<Particle, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle, false, true>::ParticleVecType *additionalVectors) const override {
    if (_isValid != ValidityState::invalid) {
      if (not particlesToAddEmpty()) {
        autopas::utils::ExceptionHandler::exception(
            "VerletClusterLists::getRegionIterator() const: Error: particle container is valid, but _particlesToAdd "
            "isn't empty!");
      }
      // If the particles are sorted into the towers, we can simply use the iteration over towers.
      return ContainerIterator<Particle, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
    } else {
      // if the particles are not sorted into the towers, we have to also iterate over _particlesToAdd.
      // store all pointers in a temporary which is passed to the ParticleIterator constructor.
      typename ContainerIterator<Particle, false, true>::ParticleVecType additionalVectorsToPass;
      if (not additionalVectorsEmpty<false, true>(additionalVectors)) {
        additionalVectorsToPass.reserve(_particlesToAdd.size() + additionalVectors->size());
        additionalVectorsToPass.insert(additionalVectorsToPass.end(), additionalVectors->begin(),
                                       additionalVectors->end());
      } else {
        additionalVectorsToPass.reserve(_particlesToAdd.size());
      }
      for (auto &vec : _particlesToAdd) {
        additionalVectorsToPass.push_back(&vec);
      }
      return ContainerIterator<Particle, false, true>(*this, behavior, &additionalVectorsToPass, lowerCorner,
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
    for (size_t i = 0; i < _towers.size(); ++i) {
      auto &tower = _towers[i];
      const auto [towerLowCorner, towerHighCorner] = getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
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
    for (size_t i = 0; i < _towers.size(); ++i) {
      auto &tower = _towers[i];
      const auto [towerLowCorner, towerHighCorner] = getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
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
    for (size_t i = 0; i < _towers.size(); ++i) {
      auto &tower = _towers[i];
      const auto [towerLowCorner, towerHighCorner] = getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
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
    for (size_t i = 0; i < _towers.size(); ++i) {
      auto &tower = _towers[i];
      const auto [towerLowCorner, towerHighCorner] = getTowerBoundingBox(i);
      // particles can move over cell borders. Calculate the volume this cell's particles can be.
      const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
      const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
      if (utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, lowerCorner, higherCorner)) {
        tower.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
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

  [[nodiscard]] unsigned long getNumberOfParticles() const override {
    size_t sum = std::accumulate(_towers.begin(), _towers.end(), 0,
                                 [](size_t acc, const auto &tower) { return acc + tower.getNumActualParticles(); });
    sum = std::accumulate(_particlesToAdd.begin(), _particlesToAdd.end(), sum,
                          [](size_t acc, const auto &buffer) { return acc + buffer.size(); });
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
   * Calculates the low and high corner of a tower given by its index.
   *
   * @note If towers are not built yet the corners of the full container are returned.
   *
   * @param index1D The tower's index in _towers.
   * @return tuple<lowCorner, highCorner>
   */
  [[nodiscard]] std::tuple<std::array<double, 3>, std::array<double, 3>> getTowerBoundingBox(size_t index1D) const {
    // case: towers are not built yet.
    if (_towersPerDim[0] == 0) {
      return {_boxMin, _boxMax};
    }
    return getTowerBoundingBox(towerIndex1DTo2D(index1D));
  }

  /**
   * Calculates the low and high corner of a tower given by its 2D grid index.
   *
   * @note If towers are not built yet the corners of the full container are returned.
   *
   * @param index2D The tower's 2D index in the grid.
   * @return tuple<lowCorner, highCorner>
   */
  [[nodiscard]] std::tuple<std::array<double, 3>, std::array<double, 3>> getTowerBoundingBox(
      const std::array<size_t, 2> &index2D) const {
    // case: towers are not built yet.
    if (_towersPerDim[0] == 0) {
      return {_boxMin, _boxMax};
    }
    const std::array<double, 3> boxMin{
        _haloBoxMin[0] + _towerSideLength * static_cast<double>(index2D[0]),
        _haloBoxMin[1] + _towerSideLength * static_cast<double>(index2D[1]),
        _haloBoxMin[2],
    };
    const std::array<double, 3> boxMax{
        boxMin[0] + _towerSideLength,
        boxMin[1] + _towerSideLength,
        _haloBoxMax[2],
    };
    return {boxMin, boxMax};
  }

  /**
   * Return the 2D index of the tower at a given position in the global domain.
   * @param pos
   * @return array<X, Y>
   */
  [[nodiscard]] std::array<size_t, 2> getTowerIndex2DAtPosition(const std::array<double, 3> &pos) const {
    // special case: Towers are not yet built, then everything is in this one tower.
    if (_towers.size() == 1) {
      return {0, 0};
    }

    // FIXME: this is a copy of VerletClusterListsRebuilder::getTowerCoordinates(). See Issue 716
    // https://github.com/AutoPas/AutoPas/issues/716
    std::array<size_t, 2> towerIndex2D{};

    for (int dim = 0; dim < 2; dim++) {
      const auto towerDimIndex =
          static_cast<long int>(floor((pos[dim] - _haloBoxMin[dim]) * _towerSideLengthReciprocal));
      const auto towerDimIndexNonNegative = static_cast<size_t>(std::max(towerDimIndex, 0l));
      const auto towerDimIndexNonLargerValue = std::min(towerDimIndexNonNegative, _towersPerDim[dim] - 1);
      towerIndex2D[dim] = towerDimIndexNonLargerValue;
      /// @todo this is a sanity check to prevent doubling of particles, but could be done better! e.g. by border and
      // flag manager
      if (pos[dim] >= _boxMax[dim]) {
        towerIndex2D[dim] = _towersPerDim[dim] - 1;
      } else if (pos[dim] < _boxMin[dim]) {
        towerIndex2D[dim] = 0;
      }
    }

    return towerIndex2D;
  }

  /**
   * Return the 1D index of the tower at a given position
   * @param pos
   * @return
   */
  [[nodiscard]] size_t getTowerIndex1DAtPosition(const std::array<double, 3> &pos) const {
    const auto [x, y] = getTowerIndex2DAtPosition(pos);
    return towerIndex2DTo1D(x, y);
  }

  /**
   * Return a reference to the tower at a given position in the simulation coordinate system (e.g. particle position).
   * @param pos
   * @return
   */
  internal::ClusterTower<Particle> &getTowerAtPosition(const std::array<double, 3> &pos) {
    const auto [x, y] = getTowerIndex2DAtPosition(pos);
    return getTowerByIndex(x, y);
  }

  /**
   * Returns a reference to the tower for the given tower grid coordinates.
   * @param x The x-th tower in x direction.
   * @param y The y-th tower in y direction.
   * @return a reference to the tower for the given tower grid coordinates.
   */
  internal::ClusterTower<Particle> &getTowerByIndex(const size_t x, const size_t y) {
    return _towers[towerIndex2DTo1D(x, y)];
  }

  /**
   * Returns the 1D index for the given tower grid coordinates of a tower.
   *
   * @param x The x-coordinate of the tower.
   * @param y The y-coordinate of the tower.
   * @param towersPerDim The number of towers in each dimension.
   * @return the 1D index for the given tower grid coordinates of a tower.
   */
  static size_t towerIndex2DTo1D(const size_t x, const size_t y, const std::array<size_t, 2> &towersPerDim) {
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
   * Returns the 2D index for the given 1D index of a tower.
   *
   * @param index
   * @return the 2D index for the given 1D index of a tower.
   */
  [[nodiscard]] std::array<size_t, 2> towerIndex1DTo2D(size_t index) const {
    if (_towersPerDim[0] == 0) {
      return {0, 0};
    } else {
      return {index % _towersPerDim[0], index / _towersPerDim[0]};
    }
  }

  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override { return _boxMax; }

  void setBoxMax(const std::array<double, 3> &boxMax) override { _boxMax = boxMax; }

  /**
   * Get the upper corner of the halo box.
   * @return the upper corner of the halo box.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMax() const { return _haloBoxMax; }

  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override { return _boxMin; }

  void setBoxMin(const std::array<double, 3> &boxMin) override { _boxMin = boxMin; }

  /**
   * Get the lower corner of the halo box.
   * @return the lower corner of the halo box.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMin() const { return _haloBoxMin; }

  [[nodiscard]] double getCutoff() const override { return _cutoff; }

  void setCutoff(double cutoff) override { _cutoff = cutoff; }

  [[nodiscard]] double getVerletSkin() const override { return _skinPerTimestep * _rebuildFrequency; }

  /**
   * Set the verlet skin length per timestep for the container.
   * @param skinPerTimestep
   */
  void setSkinPerTimestep(double skinPerTimestep) { _skinPerTimestep = skinPerTimestep; }

  /**
   * Get the rebuild Frequency value for the container.
   * @return rebuildFrequency
   */
  [[nodiscard]] unsigned int getRebuildFrequency() { return _rebuildFrequency; }

  /**
   * Set the rebuild Frequency value for the container.
   * @param rebuildFrequency
   */
  void setRebuildFrequency(unsigned int rebuildFrequency) { _rebuildFrequency = rebuildFrequency; }

  [[nodiscard]] double getInteractionLength() const override { return _cutoff + _skinPerTimestep * _rebuildFrequency; }

  void deleteAllParticles() override {
    _isValid.store(ValidityState::invalid, std::memory_order::memory_order_relaxed);
    std::for_each(_particlesToAdd.begin(), _particlesToAdd.end(), [](auto &buffer) { buffer.clear(); });
    std::for_each(_towers.begin(), _towers.end(), [](auto &tower) { tower.clear(); });
  }

 protected:
  /**
   * Rebuild the towers and the clusters.
   * This function sets the container structure to valid.
   */
  void rebuildTowersAndClusters() {
    // collect all particles to add from across the thread buffers
    typename decltype(_particlesToAdd)::value_type particlesToAdd;
    size_t numParticlesToAdd = std::accumulate(_particlesToAdd.begin(), _particlesToAdd.end(), 0,
                                               [](size_t acc, const auto &buffer) { return acc + buffer.size(); });
    particlesToAdd.reserve(numParticlesToAdd);
    std::for_each(_particlesToAdd.begin(), _particlesToAdd.end(), [&](auto &particlesBuffer) {
      particlesToAdd.insert(particlesToAdd.end(), particlesBuffer.begin(), particlesBuffer.end());
      particlesBuffer.clear();
    });

    _builder =
        std::make_unique<internal::VerletClusterListsRebuilder<Particle>>(*this, _towers, particlesToAdd, _clusterSize);

    std::tie(_towerSideLength, _numTowersPerInteractionLength, _towersPerDim, _numClusters) =
        _builder->rebuildTowersAndClusters();

    _towerSideLengthReciprocal = 1 / _towerSideLength;
    _isValid.store(ValidityState::cellsValidListsInvalid, std::memory_order::memory_order_relaxed);
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
        auto &tower = getTowerByIndex(x, y);
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
        auto &tower = getTowerByIndex(x, y);

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
   * Given a pair of cell-/particleIndex and iterator restrictions either returns the next indices that match these
   * restrictions or indices that are out of bounds (e.g. cellIndex >= cells.size())
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @param endCellIndex
   * @return tuple<cellIndex, particleIndex>
   */
  template <bool regionIter>
  [[nodiscard]] std::tuple<size_t, size_t> advanceIteratorIndices(size_t cellIndex, size_t particleIndex,
                                                                  IteratorBehavior iteratorBehavior,
                                                                  const std::array<double, 3> &boxMin,
                                                                  const std::array<double, 3> &boxMax,
                                                                  size_t endCellIndex) const {
    // Finding the indices for the next particle
    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();

    // helper function to determine if the cell can even contain particles of interest to the iterator
    auto towerIsRelevant = [&]() -> bool {
      // special case: Towers are not yet built, then the tower acting as buffer is always relevant.
      if (_towers.size() == 1) {
        return true;
      }
      bool isRelevant = true;
      if constexpr (regionIter) {
        // is the cell in the region?
        const auto [towerLowCorner, towerHighCorner] = getTowerBoundingBox(cellIndex);
        // particles can move over cell borders. Calculate the volume this cell's particles can be.
        const auto towerLowCornerSkin = utils::ArrayMath::subScalar(towerLowCorner, this->getVerletSkin() * 0.5);
        const auto towerHighCornerSkin = utils::ArrayMath::addScalar(towerHighCorner, this->getVerletSkin() * 0.5);
        isRelevant = utils::boxesOverlap(towerLowCornerSkin, towerHighCornerSkin, boxMin, boxMax);
      }
      return isRelevant;
    };

    do {
      // advance to the next particle
      ++particleIndex;
      // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.

      // If cell has wrong type, or there are no more particles in this cell jump to the next
      while (not towerIsRelevant() or particleIndex >= this->_towers[cellIndex].numParticles()) {
        cellIndex += stride;
        particleIndex = 0;

        // If we notice that there is nothing else to look at set invalid values, so we get a nullptr next time and
        // break.
        if (cellIndex > endCellIndex) {
          return {std::numeric_limits<size_t>::max(), particleIndex};
        }
      }
    } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
        this->_towers[cellIndex][particleIndex], iteratorBehavior, boxMin, boxMax));

    // the indices returned at this point should always be valid
    return {cellIndex, particleIndex};
  }

  /**
   * load estimation algorithm for balanced traversals.
   */
  autopas::LoadEstimatorOption _loadEstimator;

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
   * Side length of xy-grid cells.
   */
  double _towerSideLength{0.};

  /**
   * 1/side length of xy-grid cells.
   */
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
   * Outer vector is for Thread buffer to allow parallel particle insertion.
   */
  std::vector<std::vector<Particle>> _particlesToAdd;

  /**
   * Checks if there are particles in at least one thread buffer of _particlesToAdd.
   * @return true iff all thread buffer are empty.
   */
  [[nodiscard]] bool particlesToAddEmpty() const {
    for (auto &threadBuffer : _particlesToAdd) {
      if (not threadBuffer.empty()) {
        return false;
      }
    }
    return true;
  }

  /**
   * Helper function for begin() and getRegionIterator() to check if the passed additionalVectors are empty
   *
   * @tparam modifiable
   * @tparam regionIter
   * @param additionalVectors
   * @return true iff all passed vectors are empty
   * @return false iff at least one buffer from the passed vectors is not empty
   */
  template <bool modifiable, bool regionIter>
  bool additionalVectorsEmpty(
      typename ContainerIterator<Particle, modifiable, regionIter>::ParticleVecType *additionalVectors) const {
    if (additionalVectors) {
      for (const auto &buffer : *additionalVectors) {
        if (not buffer->empty()) {
          return false;
        }
      }
    }
    return true;
  }

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
   * SkinPerTimestep.
   */
  double _skinPerTimestep{};
  /**
   * rebuidFrequency.
   */
  unsigned int _rebuildFrequency{};
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
};

}  // namespace autopas
