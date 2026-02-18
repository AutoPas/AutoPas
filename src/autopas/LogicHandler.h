/**
 * @file LogicHandler.h
 * @author seckler
 * @date 31.05.19
 */

#pragma once
#include <atomic>
#include <limits>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <vector>

#include "autopas/LogicHandlerInfo.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/remainder/RemainderPairwiseInteractionHandler.h"
#include "autopas/remainder/RemainderTriwiseInteractionHandler.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/ContainerSelectorInfo.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/NumParticlesEstimator.h"
#include "autopas/utils/StaticContainerSelector.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/logging/FLOPLogger.h"
#include "autopas/utils/logging/IterationLogger.h"
#include "autopas/utils/logging/IterationMeasurements.h"
#include "autopas/utils/logging/LiveInfoLogger.h"
#include "autopas/utils/logging/Logger.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas {

/**
 * The LogicHandler takes care of the containers s.t. they are all in the same valid state.
 * This is mainly done by incorporating a global container rebuild frequency, which defines when containers and their
 * neighbor lists will be rebuilt.
 */
template <typename Particle_T>
class LogicHandler {
 public:
  /**
   * Constructor of the LogicHandler.
   * @param autotuners Unordered map with interaction types and respective autotuner instances.
   * @param logicHandlerInfo
   * @param rebuildFrequency
   * @param outputSuffix
   */
  LogicHandler(std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> &autotuners,
               const LogicHandlerInfo &logicHandlerInfo, unsigned int rebuildFrequency, const std::string &outputSuffix)
      : _autoTunerRefs(autotuners),
        _logicHandlerInfo(logicHandlerInfo),
        _neighborListRebuildFrequency{rebuildFrequency},
        _particleBuffer(autopas_get_max_threads()),
        _haloParticleBuffer(autopas_get_max_threads()),
        _remainderPairwiseInteractionHandler(_spatialLocks),
        _remainderTriwiseInteractionHandler(_spatialLocks),
        _verletClusterSize(logicHandlerInfo.verletClusterSize),
        _sortingThreshold(logicHandlerInfo.sortingThreshold),
        _iterationLogger(outputSuffix, std::any_of(autotuners.begin(), autotuners.end(),
                                                   [](const auto &tuner) { return tuner.second->canMeasureEnergy(); })),
        _flopLogger(outputSuffix),
        _liveInfoLogger(outputSuffix) {
    using namespace autopas::utils::ArrayMath::literals;
    // Initialize AutoPas with tuners for given interaction types
    for (const auto &[interactionType, tuner] : autotuners) {
      _interactionTypes.insert(interactionType);

      const auto configuration = tuner->getCurrentConfig();
      // initialize the container and make sure it is valid
      _currentContainerSelectorInfo = ContainerSelectorInfo{_logicHandlerInfo.boxMin,
                                                            _logicHandlerInfo.boxMax,
                                                            _logicHandlerInfo.cutoff,
                                                            configuration.cellSizeFactor,
                                                            _logicHandlerInfo.verletSkin,
                                                            _verletClusterSize,
                                                            _sortingThreshold,
                                                            configuration.loadEstimator};
      _currentContainer =
          ContainerSelector<Particle_T>::generateContainer(configuration.container, _currentContainerSelectorInfo);
      checkMinimalSize();
    }

    // initialize locks needed for remainder traversal
    const auto interactionLength = logicHandlerInfo.cutoff + logicHandlerInfo.verletSkin;
    const auto interactionLengthInv = 1. / interactionLength;
    const auto boxLengthWithHalo = logicHandlerInfo.boxMax - logicHandlerInfo.boxMin + (2 * interactionLength);
    initSpatialLocks(boxLengthWithHalo, interactionLengthInv);
  }

  /**
   * Returns a non-const reference to the currently selected particle container.
   * @return Non-const reference to the container.
   */
  ParticleContainerInterface<Particle_T> &getContainer() { return *_currentContainer; }

  /**
   * Collects leaving particles from buffer and potentially inserts owned particles to the container.
   * @param insertOwnedParticlesToContainer Decides whether to insert owned particles to the container.
   * @return Leaving particles.
   */
  [[nodiscard]] std::vector<Particle_T> collectLeavingParticlesFromBuffer(bool insertOwnedParticlesToContainer) {
    const auto &boxMin = _currentContainer->getBoxMin();
    const auto &boxMax = _currentContainer->getBoxMax();
    std::vector<Particle_T> leavingBufferParticles{};
    for (auto &cell : _particleBuffer) {
      auto &buffer = cell._particles;
      if (insertOwnedParticlesToContainer) {
        // Can't be const because we potentially modify ownership before re-adding
        for (auto &p : buffer) {
          if (p.isDummy()) {
            continue;
          }
          if (utils::inBox(p.getR(), boxMin, boxMax)) {
            p.setOwnershipState(OwnershipState::owned);
            _currentContainer->addParticle(p);
          } else {
            leavingBufferParticles.push_back(p);
          }
        }
        buffer.clear();
      } else {
        for (auto iter = buffer.begin(); iter < buffer.end();) {
          auto &p = *iter;

          auto fastRemoveP = [&]() {
            // Fast remove of particle, i.e., swap with last entry && pop.
            std::swap(p, buffer.back());
            buffer.pop_back();
            // Do not increment the iter afterward!
          };
          if (p.isDummy()) {
            // We remove dummies!
            fastRemoveP();
            // In case we swapped a dummy here, don't increment the iterator and do another iteration to check again.
            continue;
          }
          // if p was a dummy a new particle might now be at the memory location of p so we need to check that.
          // We also just might have deleted the last particle in the buffer in that case the inBox check is meaningless
          if (not buffer.empty() and utils::notInBox(p.getR(), boxMin, boxMax)) {
            leavingBufferParticles.push_back(p);
            fastRemoveP();
          } else {
            ++iter;
          }
        }
      }
    }
    return leavingBufferParticles;
  }

  /**
   * @copydoc AutoPas::updateContainer()
   */
  [[nodiscard]] std::vector<Particle_T> updateContainer() {
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    this->checkNeighborListsInvalidDoDynamicRebuild();
#endif
    bool doDataStructureUpdate = not neighborListsAreValid();

    if (_functorCalls > 0) {
      // Bump iteration counters for all autotuners
      for (const auto &[interactionType, autoTuner] : _autoTunerRefs) {
        const bool needsToWait = checkTuningStates(interactionType);
        // Called before bumpIterationCounters as it would return false after that.
        if (autoTuner->inLastTuningIteration()) {
          _iterationAtEndOfLastTuningPhase = _iteration;
        }
        autoTuner->bumpIterationCounters(needsToWait);
      }

      // We will do a rebuild in this timestep
      if (not _neighborListsAreValid.load(std::memory_order_relaxed)) {
        _stepsSinceLastListRebuild = 0;
      }
      ++_stepsSinceLastListRebuild;
      _currentContainer->setStepsSinceLastRebuild(_stepsSinceLastListRebuild);
      ++_iteration;
    }

    // The next call also adds particles to the container if doDataStructureUpdate is true.
    auto leavingBufferParticles = collectLeavingParticlesFromBuffer(doDataStructureUpdate);

    AutoPasLog(TRACE, "Initiating container update.");
    auto leavingParticles = _currentContainer->updateContainer(not doDataStructureUpdate);
    leavingParticles.insert(leavingParticles.end(), leavingBufferParticles.begin(), leavingBufferParticles.end());

    // Substract the amount of leaving particles from the number of owned particles.
    _numParticlesOwned.fetch_sub(leavingParticles.size(), std::memory_order_relaxed);
    // updateContainer deletes all halo particles.
    std::for_each(_haloParticleBuffer.begin(), _haloParticleBuffer.end(), [](auto &buffer) { buffer.clear(); });
    _numParticlesHalo.store(0, std::memory_order_relaxed);
    return leavingParticles;
  }

  /**
   * Pass values to the actual container.
   * @param boxMin
   * @param boxMax
   * @return Vector of particles that are outside the box after the resize.
   */
  std::vector<Particle_T> resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
    using namespace autopas::utils::ArrayMath::literals;
    const auto &oldMin = _currentContainer->getBoxMin();
    const auto &oldMax = _currentContainer->getBoxMax();

    // if nothing changed, do nothing
    if (oldMin == boxMin and oldMax == boxMax) {
      return {};
    }

    // sanity check that new size is actually positive
    for (size_t i = 0; i < boxMin.size(); ++i) {
      if (boxMin[i] >= boxMax[i]) {
        utils::ExceptionHandler::exception(
            "New box size in dimension {} is not positive!\nboxMin[{}] = {}\nboxMax[{}] = {}", i, i, boxMin[i], i,
            boxMax[i]);
      }
    }

    // warn if domain changes too drastically
    const auto newLength = boxMax - boxMin;
    const auto oldLength = oldMax - oldMin;
    const auto relDiffLength = newLength / oldLength;
    for (size_t i = 0; i < newLength.size(); ++i) {
      // warning threshold is set arbitrary and up for change if needed
      if (relDiffLength[i] > 1.3 or relDiffLength[i] < 0.7) {
        AutoPasLog(WARN,
                   "LogicHandler.resize(): Domain size changed drastically in dimension {}! Gathered AutoTuning "
                   "information might not be applicable anymore!\n"
                   "Size old box : {}\n"
                   "Size new box : {}\n"
                   "Relative diff: {}",
                   i, utils::ArrayUtils::to_string(oldLength), utils::ArrayUtils::to_string(newLength),
                   utils::ArrayUtils::to_string(relDiffLength));
      }
    }

    // The new box size is valid, so update the current container info.
    _currentContainerSelectorInfo.boxMin = boxMin;
    _currentContainerSelectorInfo.boxMax = boxMax;

    // check all particles
    std::vector<Particle_T> particlesNowOutside;
    for (auto pIter = _currentContainer->begin(); pIter.isValid(); ++pIter) {
      // make sure only owned ones are present
      if (not pIter->isOwned()) {
        utils::ExceptionHandler::exception(
            "LogicHandler::resizeBox() encountered non owned particle. "
            "When calling resizeBox() these should be already deleted. "
            "This could be solved by calling updateContainer() before resizeBox().");
      }
      // owned particles that are now outside are removed from the container and returned
      if (not utils::inBox(pIter->getR(), boxMin, boxMax)) {
        particlesNowOutside.push_back(*pIter);
        decreaseParticleCounter(*pIter);
        internal::markParticleAsDeleted(*pIter);
      }
    }

    // Resize by generating a new container with the new box size and moving all particles to it.
    auto newContainer = ContainerSelector<Particle_T>::generateContainer(_currentContainer->getContainerType(),
                                                                         _currentContainerSelectorInfo);
    setCurrentContainer(std::move(newContainer));
    // The container might have changed sufficiently that we would need a different number of spatial locks.
    const auto boxLength = boxMax - boxMin;
    const auto interactionLengthInv = 1. / _currentContainer->getInteractionLength();
    initSpatialLocks(boxLength, interactionLengthInv);

    // Set this flag, s.t., the container is rebuilt!
    _neighborListsAreValid.store(false, std::memory_order_relaxed);

    return particlesNowOutside;
  }

  /**
   * Estimates the number of halo particles via autopas::utils::NumParticlesEstimator::estimateNumHalosUniform(), then
   * calls LogicHandler::reserve(size_t numParticles, size_t numHaloParticles).
   *
   * @param numParticles Total number of owned particles.
   */
  void reserve(size_t numParticles) {
    const auto numParticlesHaloEstimate = autopas::utils::NumParticlesEstimator::estimateNumHalosUniform(
        numParticles, _currentContainer->getBoxMin(), _currentContainer->getBoxMax(),
        _currentContainer->getInteractionLength());
    reserve(numParticles, numParticlesHaloEstimate);
  }

  /**
   * Reserves space in the particle buffers and the container.
   *
   * @param numParticles Total number of owned particles.
   * @param numHaloParticles Total number of halo particles.
   */
  void reserve(size_t numParticles, size_t numHaloParticles) {
    const auto numHaloParticlesPerBuffer = numHaloParticles / _haloParticleBuffer.size();
    for (auto &buffer : _haloParticleBuffer) {
      buffer.reserve(numHaloParticlesPerBuffer);
    }
    // there is currently no good heuristic for this buffer, so reuse the one for halos.
    for (auto &buffer : _particleBuffer) {
      buffer.reserve(numHaloParticlesPerBuffer);
    }

    // reserve is called for the container only in the rebuild iterations.
    // during non-rebuild iterations, particles are not added in the container but in buffer.
    if (not _neighborListsAreValid.load(std::memory_order_relaxed)) {
      _currentContainer->reserve(numParticles, numHaloParticles);
    }
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(const Particle_T &p) {
    // first check that the particle actually belongs in the container
    const auto &boxMin = _currentContainer->getBoxMin();
    const auto &boxMax = _currentContainer->getBoxMax();
    if (utils::notInBox(p.getR(), boxMin, boxMax)) {
      autopas::utils::ExceptionHandler::exception(
          "LogicHandler: Trying to add a particle that is not in the bounding box.\n"
          "Box Min {}\n"
          "Box Max {}\n"
          "{}",
          boxMin, boxMax, p.toString());
    }
    Particle_T particleCopy = p;
    particleCopy.setOwnershipState(OwnershipState::owned);
    if (not _neighborListsAreValid.load(std::memory_order_relaxed)) {
      // Container has to (about to) be invalid to be able to add Particles!
      _currentContainer->template addParticle<false>(particleCopy);
    } else {
      // If the container is valid, we add it to the particle buffer.
      _particleBuffer[autopas_get_thread_num()].addParticle(particleCopy);
    }
    _numParticlesOwned.fetch_add(1, std::memory_order_relaxed);
  }

  /**
   * @copydoc AutoPas::addHaloParticle()
   */
  void addHaloParticle(const Particle_T &haloParticle) {
    const auto &boxMin = _currentContainer->getBoxMin();
    const auto &boxMax = _currentContainer->getBoxMax();
    Particle_T haloParticleCopy = haloParticle;
    if (utils::inBox(haloParticleCopy.getR(), boxMin, boxMax)) {
      utils::ExceptionHandler::exception(
          "LogicHandler: Trying to add a halo particle that is not outside the box of the container.\n"
          "Box Min {}\n"
          "Box Max {}\n"
          "{}",
          utils::ArrayUtils::to_string(boxMin), utils::ArrayUtils::to_string(boxMax), haloParticleCopy.toString());
    }
    haloParticleCopy.setOwnershipState(OwnershipState::halo);
    if (not _neighborListsAreValid.load(std::memory_order_relaxed)) {
      // If the neighbor lists are not valid, we can add the particle.
      _currentContainer->template addHaloParticle</* checkInBox */ false>(haloParticleCopy);
    } else {
      // Check if we can update an existing halo(dummy) particle.
      bool updated = _currentContainer->updateHaloParticle(haloParticleCopy);
      if (not updated) {
        // If we couldn't find an existing particle, add it to the halo particle buffer.
        _haloParticleBuffer[autopas_get_thread_num()].addParticle(haloParticleCopy);
      }
    }
    _numParticlesHalo.fetch_add(1, std::memory_order_relaxed);
  }

  /**
   * @copydoc AutoPas::deleteAllParticles()
   */
  void deleteAllParticles() {
    _neighborListsAreValid.store(false, std::memory_order_relaxed);
    _currentContainer->deleteAllParticles();
    std::for_each(_particleBuffer.begin(), _particleBuffer.end(), [](auto &buffer) { buffer.clear(); });
    std::for_each(_haloParticleBuffer.begin(), _haloParticleBuffer.end(), [](auto &buffer) { buffer.clear(); });
    // all particles are gone -> reset counters.
    _numParticlesOwned.store(0, std::memory_order_relaxed);
    _numParticlesHalo.store(0, std::memory_order_relaxed);
  }

  /**
   * Takes a particle, checks if it is in any of the particle buffers, and deletes it from them if found.
   * @param particle Particle to delete. If something was deleted this reference might point to a different particle or
   * invalid memory.
   * @return Tuple: <True iff the particle was found and deleted, True iff the reference is valid>
   */
  std::tuple<bool, bool> deleteParticleFromBuffers(Particle_T &particle) {
    // find the buffer the particle belongs to
    auto &bufferCollection = particle.isOwned() ? _particleBuffer : _haloParticleBuffer;
    for (auto &cell : bufferCollection) {
      auto &buffer = cell._particles;
      // if the address of the particle is between start and end of the buffer it is in this buffer
      if (not buffer.empty() and &(buffer.front()) <= &particle and &particle <= &(buffer.back())) {
        const bool isRearParticle = &particle == &buffer.back();
        // swap-delete
        particle = buffer.back();
        buffer.pop_back();
        return {true, not isRearParticle};
      }
    }
    return {false, true};
  }

  /**
   * Decrease the correct internal particle counters.
   * This function should always be called if individual particles are deleted.
   * @param particle reference to particles that should be deleted
   */
  void decreaseParticleCounter(Particle_T &particle) {
    if (particle.isOwned()) {
      _numParticlesOwned.fetch_sub(1, std::memory_order_relaxed);
    } else {
      _numParticlesHalo.fetch_sub(1, std::memory_order_relaxed);
    }
  }

  /**
   * This function covers the full pipeline of all mechanics happening during the computation of particle interactions.
   * This includes:
   * - selecting a configuration
   *   - gather live info, homogeneity, and max density
   *   - get next config (tuning)
   *   - check applicability
   *   - instantiation of traversal and container
   * - triggering iteration and tuning result logger
   * - computing the interactions
   *   - init and end traversal
   *   - remainder traversal
   *   - measurements
   * - pass measurements to tuner
   *
   * @tparam Functor
   * @param functor
   * @param interactionType
   * @return True if this was a tuning iteration.
   */
  template <class Functor>
  bool computeInteractionsPipeline(Functor *functor, const InteractionTypeOption &interactionType);

  /**
   * Create the additional vectors vector for a given iterator behavior.
   * @tparam Iterator
   * @param behavior
   * @return Vector of pointers to buffer vectors.
   */
  template <class Iterator>
  typename Iterator::ParticleVecType gatherAdditionalVectors(IteratorBehavior behavior) {
    typename Iterator::ParticleVecType additionalVectors;
    if (not(behavior & IteratorBehavior::containerOnly)) {
      additionalVectors.reserve(static_cast<bool>(behavior & IteratorBehavior::owned) * _particleBuffer.size() +
                                static_cast<bool>(behavior & IteratorBehavior::halo) * _haloParticleBuffer.size());
      if (behavior & IteratorBehavior::owned) {
        for (auto &buffer : _particleBuffer) {
          // Don't insert empty buffers. This also means that we won't pick up particles added during iterating if they
          // go to the buffers. But since we wouldn't pick them up if they go into the container to a cell that the
          // iterators already passed this is unsupported anyways.
          if (not buffer.isEmpty()) {
            additionalVectors.push_back(&(buffer._particles));
          }
        }
      }
      if (behavior & IteratorBehavior::halo) {
        for (auto &buffer : _haloParticleBuffer) {
          if (not buffer.isEmpty()) {
            additionalVectors.push_back(&(buffer._particles));
          }
        }
      }
    }
    return additionalVectors;
  }

  /**
   * @copydoc AutoPas::begin()
   */
  ContainerIterator<Particle_T, true, false> begin(IteratorBehavior behavior) {
    auto additionalVectors = gatherAdditionalVectors<ContainerIterator<Particle_T, true, false>>(behavior);
    return _currentContainer->begin(behavior, std::ref(additionalVectors));
  }

  /**
   * @copydoc AutoPas::begin()
   */
  ContainerIterator<Particle_T, false, false> begin(IteratorBehavior behavior) const {
    auto additionalVectors =
        const_cast<LogicHandler *>(this)->gatherAdditionalVectors<ContainerIterator<Particle_T, false, false>>(
            behavior);
    return _currentContainer->begin(behavior, std::ref(additionalVectors));
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  ContainerIterator<Particle_T, true, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                              const std::array<double, 3> &higherCorner,
                                                              IteratorBehavior behavior) {
    // sanity check: Most of our stuff depends on `inBox`, which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner[d] > higherCorner[d]) {
        utils::ExceptionHandler::exception(
            "Requesting region Iterator where the upper corner is lower than the lower corner!\n"
            "Lower corner: {}\n"
            "Upper corner: {}",
            lowerCorner, higherCorner);
      }
    }

    auto additionalVectors = gatherAdditionalVectors<ContainerIterator<Particle_T, true, true>>(behavior);
    return _currentContainer->getRegionIterator(lowerCorner, higherCorner, behavior, std::ref(additionalVectors));
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  ContainerIterator<Particle_T, false, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                               const std::array<double, 3> &higherCorner,
                                                               IteratorBehavior behavior) const {
    // sanity check: Most of our stuff depends on `inBox`, which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner[d] > higherCorner[d]) {
        utils::ExceptionHandler::exception(
            "Requesting region Iterator where the upper corner is lower than the lower corner!\n"
            "Lower corner: {}\n"
            "Upper corner: {}",
            lowerCorner, higherCorner);
      }
    }

    auto additionalVectors =
        const_cast<LogicHandler *>(this)->gatherAdditionalVectors<ContainerIterator<Particle_T, false, true>>(behavior);
    return std::as_const(_currentContainer)
        ->getRegionIterator(lowerCorner, higherCorner, behavior, std::ref(additionalVectors));
  }

  /**
   * Get the number of owned particles.
   * @return
   */
  [[nodiscard]] unsigned long getNumberOfParticlesOwned() const { return _numParticlesOwned; }

  /**
   * Get the number of halo particles.
   * @return
   */
  [[nodiscard]] unsigned long getNumberOfParticlesHalo() const { return _numParticlesHalo; }

  /**
   * Check if other autotuners for any other interaction types are still in a tuning phase.
   * @param interactionType
   * @return bool whether other tuners are still tuning.
   */
  bool checkTuningStates(const InteractionTypeOption &interactionType) {
    // Goes over all pairs in _autoTunerRefs and returns true as soon as one is `inTuningPhase()`.
    // The tuner associated with the given interaction type is ignored.
    return std::any_of(std::begin(_autoTunerRefs), std::end(_autoTunerRefs), [&](const auto &entry) {
      return not(entry.first == interactionType) and entry.second->inTuningPhase();
    });
  }

  /**
   * Checks if the given configuration can be used with the given functor and the current state of the simulation.
   * For this, the container and traversal need to be instantiated, hence if the configuration is applicable, it sets
   * the current container and returns the traversal.
   *
   * @tparam Functor
   * @param config
   * @param functor
   * @return tuple<Traversal, rejectIndefinitely> Traversal is a nullptr if the configuration is not applicable
   * The bool rejectIndefinitely indicates if the configuration can be completely removed from the search space because
   * it will never be applicable.
   */
  template <class Functor>
  [[nodiscard]] std::tuple<std::unique_ptr<TraversalInterface>, bool> isConfigurationApplicable(
      const Configuration &config, Functor &functor);

  /**
   * Directly exchange the internal particle and halo buffers with the given vectors and update particle counters.
   *
   * @note This function is for testing purposes only!
   * @warning This function only sets as many buffers as are given to it. E.g. if particleBuffers.size() == 3 but there
   * LogicHandler has 8 the last five buffers will not be touched.
   *
   * @param particleBuffers
   * @param haloParticleBuffers
   */
  void setParticleBuffers(const std::vector<FullParticleCell<Particle_T>> &particleBuffers,
                          const std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers);

  /**
   * Getter for the particle buffers.
   *
   * @note Intended for tests only.
   *
   * @return tuple of const references to the internal buffers.
   */
  std::tuple<const std::vector<FullParticleCell<Particle_T>> &, const std::vector<FullParticleCell<Particle_T>> &>
  getParticleBuffers() const;

  /**
   * Getter for the mean rebuild frequency.
   * Helpful for determining the frequency for the dynamic containers
   * This function is only used for dynamic containers currently but returns user defined rebuild frequency for static
   * case for safety.
   * @param considerOnlyLastNonTuningPhase Bool to determine if mean rebuild frequency is to be calculated over entire
   * iterations or only during the last non-tuning phase.
   * The mean rebuild frequency over the non-tuning phase is required by the autoTuner for weighting the rebuild and
   * non-rebuild samples.
   * @return value of the mean frequency as double
   */
  [[nodiscard]] double getMeanRebuildFrequency(bool considerOnlyLastNonTuningPhase = false) const {
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    const auto numRebuilds = considerOnlyLastNonTuningPhase ? _numRebuildsInNonTuningPhase : _numRebuilds;
    // The total number of iterations is iteration + 1
    const auto iterationCount =
        considerOnlyLastNonTuningPhase ? _iteration - _iterationAtEndOfLastTuningPhase : _iteration + 1;
    if (numRebuilds == 0) {
      return static_cast<double>(_neighborListRebuildFrequency);
    } else {
      return static_cast<double>(iterationCount) / numRebuilds;
    }
#else
    return static_cast<double>(_neighborListRebuildFrequency);
#endif
  }

  /**
   * Estimates the rebuild frequency based on the current maximum velocity in the container
   * Using the formula rf = skin/deltaT/vmax/2
   * @param skin is the skin length used in the simulation
   * @param deltaT is the time step
   * @return estimate of the current rebuild frequency
   */
  double getVelocityMethodRFEstimate(const double skin, const double deltaT) const {
    using autopas::utils::ArrayMath::dot;
    // Initialize the maximum velocity to zero
    double maxVelocity = 0;
    // Iterate over the owned particles in container to determine maximum velocity
    AUTOPAS_OPENMP(parallel reduction(max : maxVelocity))
    for (auto iter = this->begin(IteratorBehavior::owned | IteratorBehavior::containerOnly); iter.isValid(); ++iter) {
      std::array<double, 3> tempVel = iter->getV();
      double tempVelAbs = sqrt(dot(tempVel, tempVel));
      maxVelocity = std::max(tempVelAbs, maxVelocity);
    }
    // return the rebuild frequency estimate
    return skin / maxVelocity / deltaT / 2;
  }
  /**
   * getter function for _neighborListInvalidDoDynamicRebuild
   * @return bool stored in _neighborListInvalidDoDynamicRebuild
   */
  bool getNeighborListsInvalidDoDynamicRebuild();

  /**
   * Checks if any particle has moved more than skin/2.
   * updates bool: _neighborListInvalidDoDynamicRebuild
   */
  void checkNeighborListsInvalidDoDynamicRebuild();

  /**
   * Checks if any particle has moved more than skin/2.
   * resets bool: _neighborListInvalidDoDynamicRebuild to false
   */
  void resetNeighborListsInvalidDoDynamicRebuild();

  /**
   * Checks if in the next iteration the neighbor lists have to be rebuilt.
   *
   * This can be the case either because we hit the rebuild frequency or the dynamic rebuild criteria or because the
   * auto tuner tests a new configuration.
   *
   * @return True iff the neighbor lists will not be rebuild.
   */
  bool neighborListsAreValid();

 private:
  /**
   * Initialize or update the spatial locks used during the remainder traversal.
   * If the locks are already initialized but the container size changed, surplus locks will
   * be deleted, new locks are allocated and locks that are still necessary are reused.
   *
   * @note Should the generated number of locks exceed 1e6, the number of locks is reduced so that
   * the number of locks per dimensions is proportional to the domain side lengths and smaller than the limit.
   *
   * @param boxLength Size per dimension for the box that should be covered by locks (should include halo).
   * @param interactionLengthInv Inverse of the side length of the virtual boxes one lock is responsible for.
   *
   */
  void initSpatialLocks(const std::array<double, 3> &boxLength, double interactionLengthInv) {
    using namespace autopas::utils::ArrayMath::literals;
    using utils::ArrayMath::ceil;
    using utils::ArrayUtils::static_cast_copy_array;

    // The maximum number of spatial locks is capped at 1e6.
    // This limit is chosen more or less arbitrary. It is big enough so that our regular MD simulations
    // fall well within it and small enough so that no memory issues arise.
    // There were no rigorous tests for an optimal number of locks.
    // Without this cap, very large domains (or tiny cutoffs) would generate an insane number of locks,
    // that could blow up the memory.
    constexpr size_t maxNumSpacialLocks{1000000};

    // One lock per interaction length or less if this would generate too many.
    const std::array<size_t, 3> locksPerDim = [&]() {
      // First naively calculate the number of locks if we simply take the desired cell length.
      // Ceil because both decisions are possible, and we are generous gods.
      const std::array<size_t, 3> locksPerDimNaive =
          static_cast_copy_array<size_t>(ceil(boxLength * interactionLengthInv));
      const auto totalLocksNaive =
          std::accumulate(locksPerDimNaive.begin(), locksPerDimNaive.end(), 1ul, std::multiplies<>());
      // If the number of locks is within the limits everything is fine and we can return.
      if (totalLocksNaive <= maxNumSpacialLocks) {
        return locksPerDimNaive;
      } else {
        // If the number of locks grows too large, calculate the locks per dimension proportionally to the side lengths.
        // Calculate side length relative to dimension 0.
        const std::array<double, 3> boxSideProportions = {
            1.,
            boxLength[0] / boxLength[1],
            boxLength[0] / boxLength[2],
        };
        // With this, calculate the number of locks the first dimension should receive.
        const auto prodProportions =
            std::accumulate(boxSideProportions.begin(), boxSideProportions.end(), 1., std::multiplies<>());
        // Needs floor, otherwise we exceed the limit.
        const auto locksInFirstDimFloat = std::floor(std::cbrt(maxNumSpacialLocks * prodProportions));
        // From this and the proportions relative to the first dimension, we can calculate the remaining number of locks
        const std::array<size_t, 3> locksPerDimLimited = {
            static_cast<size_t>(locksInFirstDimFloat),  // omitted div by 1
            static_cast<size_t>(locksInFirstDimFloat / boxSideProportions[1]),
            static_cast<size_t>(locksInFirstDimFloat / boxSideProportions[2]),
        };
        return locksPerDimLimited;
      }
    }();
    _spatialLocks.resize(locksPerDim[0]);
    for (auto &lockVecVec : _spatialLocks) {
      lockVecVec.resize(locksPerDim[1]);
      for (auto &lockVec : lockVecVec) {
        lockVec.resize(locksPerDim[2]);
        for (auto &lockPtr : lockVec) {
          if (not lockPtr) {
            lockPtr = std::make_unique<std::mutex>();
          }
        }
      }
    }
  }

  /**
   * Gathers dynamic data from the domain if necessary and retrieves the next configuration to use.
   * @tparam Functor
   * @param functor
   * @param interactionType
   * @return
   */
  template <class Functor>
  std::tuple<Configuration, std::unique_ptr<TraversalInterface>, bool> selectConfiguration(
      Functor &functor, const InteractionTypeOption &interactionType);

  /**
   * Sets the current container to the given one and transfers all particles from the old container to the new one.
   * @param newContainer The new container to set.
   */
  void setCurrentContainer(std::unique_ptr<ParticleContainerInterface<Particle_T>> newContainer);

  /**
   * Triggers the core steps of computing the particle interactions:
   *    - functor init- / end traversal
   *    - rebuilding of neighbor lists
   *    - container.computeInteractions()
   *    - remainder traversal
   *    - time and energy measurements.
   *
   * @tparam Functor
   * @param functor
   * @param traversal
   * @return Struct containing time and energy measurements. If no energy measurements were possible the respective
   * fields are filled with NaN.
   */
  template <class Functor>
  IterationMeasurements computeInteractions(Functor &functor, TraversalInterface &traversal);

  /**
   * Select the right Remainder function depending on the interaction type and newton3 setting.
   *
   * @tparam Functor
   * @param functor
   * @param newton3
   * @param useSoA Use SoA functor calls where it is feasible. This still leaves some interactions with the AoS functor
   * but the switch is useful to disable SoA calls if a functor doesn't support them at all.
   * @return
   */
  template <class Functor>
  void computeRemainderInteractions(Functor &functor, bool newton3, bool useSoA);

  /**
   * Check that the simulation box is at least of interaction length in each direction.
   *
   * Throws an exception if minimal size requirements are violated.
   */
  void checkMinimalSize() const;

  const LogicHandlerInfo _logicHandlerInfo;
  /**
   * Specifies after how many pair-wise traversals the neighbor lists (if they exist) are to be rebuild.
   */
  unsigned int _neighborListRebuildFrequency;

  /**
   * Number of particles in a VCL cluster.
   */
  unsigned int _verletClusterSize;

  /**
   * This is used to store the total number of neighbour lists rebuild.
   */
  size_t _numRebuilds{0};

  /**
   * This is used to store the total number of neighbour lists rebuilds in the non-tuning phase.
   * This is reset at the start of a new tuning phase.
   */
  size_t _numRebuildsInNonTuningPhase{0};

  /**
   * Number of particles in two cells from which sorting should be performed for traversal that use the CellFunctor
   */
  size_t _sortingThreshold;

  /**
   * Reference to the map of AutoTuners which are managed by the AutoPas main interface.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> &_autoTunerRefs;

  /**
   * The current container holding the particles.
   */
  std::unique_ptr<ParticleContainerInterface<Particle_T>> _currentContainer{nullptr};

  /**
   * The current configuration for the container.
   */
  ContainerSelectorInfo _currentContainerSelectorInfo;

  /**
   * Handles pairwise interactions of buffer particles (the remainder).
   */
  RemainderPairwiseInteractionHandler<Particle_T> _remainderPairwiseInteractionHandler;

  /**
   * Handles triwise interactions of buffer particles (the remainder).
   */
  RemainderTriwiseInteractionHandler<Particle_T> _remainderTriwiseInteractionHandler;

  /**
   * Set of interaction types AutoPas is initialized to, determined by the given AutoTuners.
   */
  std::set<InteractionTypeOption> _interactionTypes{};

  /**
   * Specifies if the neighbor list is valid.
   */
  std::atomic<bool> _neighborListsAreValid{false};

  /**
   * Steps since last rebuild
   */
  unsigned int _stepsSinceLastListRebuild{0};

  /**
   * Total number of functor calls of all interaction types.
   */
  unsigned int _functorCalls{0};

  /**
   * The current iteration number.
   */
  unsigned int _iteration{0};

  /**
   * The iteration number at the end of last tuning phase.
   */
  unsigned int _iterationAtEndOfLastTuningPhase{0};

  /**
   * Atomic tracker of the number of owned particles.
   */
  std::atomic<size_t> _numParticlesOwned{0ul};

  /**
   * Atomic tracker of the number of halo particles.
   */
  std::atomic<size_t> _numParticlesHalo{0ul};

  /**
   * Buffer to store particles that should not yet be added to the container. There is one buffer per thread.
   */
  std::vector<FullParticleCell<Particle_T>> _particleBuffer;

  /**
   * Buffer to store halo particles that should not yet be added to the container. There is one buffer per thread.
   */
  std::vector<FullParticleCell<Particle_T>> _haloParticleBuffer;

  /**
   * Locks for regions in the domain. Used for buffer <-> container interaction.
   * The number of locks depends on the size of the domain and interaction length but is capped to avoid
   * having an insane number of locks. For details see initSpacialLocks().
   */
  std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> _spatialLocks;

  /**
   * Logger for configuration used and time spent breakdown of iteratePairwise.
   */
  IterationLogger _iterationLogger;

  /**
   * Tells if dynamic rebuild is necessary currently
   * neighborListInvalidDoDynamicRebuild - true if a particle has moved more than skin/2
   */
  bool _neighborListInvalidDoDynamicRebuild{false};

  /**
   * updating position at rebuild for every particle
   */
  void updateRebuildPositions();

  /**
   * Logger for live info
   */
  LiveInfoLogger _liveInfoLogger;

  /**
   * Logger for FLOP count and hit rate.
   */
  FLOPLogger _flopLogger;
};

template <typename Particle_T>
void LogicHandler<Particle_T>::updateRebuildPositions() {
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  // The owned particles in buffer are ignored because they do not rely on the structure of the particle containers,
  // e.g. neighbour list, and these are iterated over using the region iterator. Movement of particles in buffer doesn't
  // require a rebuild of neighbor lists.
  AUTOPAS_OPENMP(parallel)
  for (auto iter = this->begin(IteratorBehavior::owned | IteratorBehavior::containerOnly); iter.isValid(); ++iter) {
    iter->resetRAtRebuild();
  }
#endif
}

template <typename Particle_T>
void LogicHandler<Particle_T>::checkMinimalSize() const {
  // check boxSize at least cutoff + skin
  for (unsigned int dim = 0; dim < 3; ++dim) {
    if (_currentContainer->getBoxMax()[dim] - _currentContainer->getBoxMin()[dim] <
        _currentContainer->getInteractionLength()) {
      utils::ExceptionHandler::exception(
          "Box (boxMin[{}]={} and boxMax[{}]={}) is too small.\nHas to be at least cutoff({}) + skin({}) = {}.", dim,
          _currentContainer->getBoxMin()[dim], dim, _currentContainer->getBoxMax()[dim], _currentContainer->getCutoff(),
          _currentContainer->getVerletSkin(), _currentContainer->getCutoff() + _currentContainer->getVerletSkin());
    }
  }
}

template <typename Particle_T>
bool LogicHandler<Particle_T>::getNeighborListsInvalidDoDynamicRebuild() {
  return _neighborListInvalidDoDynamicRebuild;
}

template <typename Particle_T>
bool LogicHandler<Particle_T>::neighborListsAreValid() {
  // Implement rebuild indicator as function, so it is only evaluated when needed.
  const auto needRebuild = [&](const InteractionTypeOption &interactionOption) {
    return _interactionTypes.count(interactionOption) != 0 and
           _autoTunerRefs[interactionOption]->willRebuildNeighborLists();
  };

  if (_stepsSinceLastListRebuild >= _neighborListRebuildFrequency
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
      or getNeighborListsInvalidDoDynamicRebuild()
#endif
      or needRebuild(InteractionTypeOption::pairwise) or needRebuild(InteractionTypeOption::triwise)) {
    _neighborListsAreValid.store(false, std::memory_order_relaxed);
  }

  return _neighborListsAreValid.load(std::memory_order_relaxed);
}

template <typename Particle_T>
void LogicHandler<Particle_T>::checkNeighborListsInvalidDoDynamicRebuild() {
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  const auto skin = getContainer().getVerletSkin();
  // (skin/2)^2
  const auto halfSkinSquare = skin * skin * 0.25;
  // The owned particles in buffer are ignored because they do not rely on the structure of the particle containers,
  // e.g. neighbour list, and these are iterated over using the region iterator. Movement of particles in buffer doesn't
  // require a rebuild of neighbor lists.
  AUTOPAS_OPENMP(parallel reduction(or : _neighborListInvalidDoDynamicRebuild))
  for (auto iter = this->begin(IteratorBehavior::owned | IteratorBehavior::containerOnly); iter.isValid(); ++iter) {
    const auto distance = iter->calculateDisplacementSinceRebuild();
    const double distanceSquare = utils::ArrayMath::dot(distance, distance);

    _neighborListInvalidDoDynamicRebuild |= distanceSquare >= halfSkinSquare;
  }
#endif
}

template <typename Particle_T>
void LogicHandler<Particle_T>::resetNeighborListsInvalidDoDynamicRebuild() {
  _neighborListInvalidDoDynamicRebuild = false;
}

template <typename Particle_T>
void LogicHandler<Particle_T>::setParticleBuffers(
    const std::vector<FullParticleCell<Particle_T>> &particleBuffers,
    const std::vector<FullParticleCell<Particle_T>> &haloParticleBuffers) {
  auto exchangeBuffer = [](const auto &newBuffers, auto &oldBuffers, auto &particleCounter) {
    // sanity check
    if (oldBuffers.size() < newBuffers.size()) {
      autopas::utils::ExceptionHandler::exception(
          "The number of new buffers ({}) is larger than number of existing buffers ({})!", newBuffers.size(),
          oldBuffers.size());
    }

    // we will clear the old buffers so subtract the particles from the counters.
    const auto numParticlesInOldBuffers =
        std::transform_reduce(oldBuffers.begin(), std::next(oldBuffers.begin(), newBuffers.size()), 0, std::plus<>(),
                              [](const auto &cell) { return cell.size(); });
    particleCounter.fetch_sub(numParticlesInOldBuffers, std::memory_order_relaxed);

    // clear the old buffers and copy the content of the new buffers over.
    size_t numParticlesInNewBuffers = 0;
    for (size_t i = 0; i < newBuffers.size(); ++i) {
      oldBuffers[i].clear();
      for (const auto &p : newBuffers[i]) {
        ++numParticlesInNewBuffers;
        oldBuffers[i].addParticle(p);
      }
    }
    // update the counters.
    particleCounter.fetch_add(numParticlesInNewBuffers, std::memory_order_relaxed);
  };

  exchangeBuffer(particleBuffers, _particleBuffer, _numParticlesOwned);
  exchangeBuffer(haloParticleBuffers, _haloParticleBuffer, _numParticlesHalo);
}

template <typename Particle_T>
std::tuple<const std::vector<FullParticleCell<Particle_T>> &, const std::vector<FullParticleCell<Particle_T>> &>
LogicHandler<Particle_T>::getParticleBuffers() const {
  return {_particleBuffer, _haloParticleBuffer};
}

template <typename Particle_T>
template <class Functor>
IterationMeasurements LogicHandler<Particle_T>::computeInteractions(Functor &functor, TraversalInterface &traversal) {
  // Helper to derive the Functor type at compile time
  constexpr auto interactionType = [] {
    if (utils::isPairwiseFunctor<Functor>()) {
      return InteractionTypeOption::pairwise;
    } else if (utils::isTriwiseFunctor<Functor>()) {
      return InteractionTypeOption::triwise;
    } else {
      utils::ExceptionHandler::exception(
          "LogicHandler::computeInteractions(): Functor is not valid. Only pairwise and triwise functors are "
          "supported. "
          "Please use a functor derived from "
          "PairwiseFunctor or TriwiseFunctor.");
    }
  }();

  auto &autoTuner = *_autoTunerRefs[interactionType];
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  if (autoTuner.inFirstTuningIteration()) {
    _numRebuildsInNonTuningPhase = 0;
  }

  // Rebuild frequency estimation should be triggered early in the tuning phase.
  // This is necessary because runtime prediction for each trial configuration
  // depends on the rebuild frequency.
  // To avoid the influence of poorly initialized velocities at the start of the simulation,
  // the rebuild frequency is estimated at iteration corresponding to the last sample of the first configuration.
  // The rebuild frequency estimated here is then reused for the remainder of the tuning phase.
  if (autoTuner.inFirstConfigurationLastSample()) {
    // Fetch the needed information for estimating the rebuild frequency from _logicHandlerInfo
    // and estimate the current rebuild frequency using the velocity method.
    double rebuildFrequencyEstimate =
        getVelocityMethodRFEstimate(_logicHandlerInfo.verletSkin, _logicHandlerInfo.deltaT);
    double userProvidedRF = static_cast<double>(_neighborListRebuildFrequency);
    // The user defined rebuild frequency is considered as the upper bound.
    // If velocity method estimate exceeds upper bound, set the rebuild frequency to the user defined value.
    // This is done because we currently use the user defined rebuild frequency as the upper bound to avoid expensive
    // buffer interactions.
    if (rebuildFrequencyEstimate > userProvidedRF) {
      autoTuner.setRebuildFrequency(userProvidedRF);
    } else {
      autoTuner.setRebuildFrequency(rebuildFrequencyEstimate);
    }
  }
#endif
  utils::Timer timerTotal;
  utils::Timer timerRebuild;
  utils::Timer timerComputeInteractions;
  utils::Timer timerComputeRemainder;

  const bool energyMeasurementsPossible = autoTuner.resetEnergy();

  timerTotal.start();
  functor.initTraversal();

  // if lists are not valid -> rebuild;
  if (not _neighborListsAreValid.load(std::memory_order_relaxed)) {
    timerRebuild.start();
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    this->updateRebuildPositions();
#endif
    _currentContainer->rebuildNeighborLists(&traversal);
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    this->resetNeighborListsInvalidDoDynamicRebuild();
    _numRebuilds++;
    if (not autoTuner.inTuningPhase()) {
      _numRebuildsInNonTuningPhase++;
    }
#endif
    timerRebuild.stop();
    _neighborListsAreValid.store(true, std::memory_order_relaxed);
  }

  timerComputeInteractions.start();
  _currentContainer->computeInteractions(&traversal);
  timerComputeInteractions.stop();

  timerComputeRemainder.start();
  const bool newton3 = autoTuner.getCurrentConfig().newton3;
  const auto dataLayout = autoTuner.getCurrentConfig().dataLayout;
  computeRemainderInteractions(functor, newton3, dataLayout);
  timerComputeRemainder.stop();

  functor.endTraversal(newton3);

  const auto [energyWatts, energyJoules, energyDeltaT, energyTotal] = autoTuner.sampleEnergy();
  timerTotal.stop();

  constexpr auto nanD = std::numeric_limits<double>::quiet_NaN();
  constexpr auto nanL = std::numeric_limits<long>::quiet_NaN();
  return {timerComputeInteractions.getTotalTime(),
          timerComputeRemainder.getTotalTime(),
          timerRebuild.getTotalTime(),
          timerTotal.getTotalTime(),
          energyMeasurementsPossible,
          energyMeasurementsPossible ? energyWatts : nanD,
          energyMeasurementsPossible ? energyJoules : nanD,
          energyMeasurementsPossible ? energyDeltaT : nanD,
          energyMeasurementsPossible ? energyTotal : nanL};
}

template <typename Particle_T>
template <class Functor>
void LogicHandler<Particle_T>::computeRemainderInteractions(Functor &functor, bool newton3, bool useSoA) {
  withStaticContainerType(*_currentContainer, [&](auto &actualContainerType) {
    if constexpr (utils::isPairwiseFunctor<Functor>()) {
      if (newton3) {
        _remainderPairwiseInteractionHandler.template computeRemainderInteractions<true>(
            &functor, actualContainerType, _particleBuffer, _haloParticleBuffer, useSoA);
      } else {
        _remainderPairwiseInteractionHandler.template computeRemainderInteractions<false>(
            &functor, actualContainerType, _particleBuffer, _haloParticleBuffer, useSoA);
      }
    } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
      if (newton3) {
        _remainderTriwiseInteractionHandler.template computeRemainderInteractions<true>(
            &functor, actualContainerType, _particleBuffer, _haloParticleBuffer);
      } else {
        _remainderTriwiseInteractionHandler.template computeRemainderInteractions<false>(
            &functor, actualContainerType, _particleBuffer, _haloParticleBuffer);
      }
    }
  });
}

template <typename Particle_T>
template <class Functor>
std::tuple<Configuration, std::unique_ptr<TraversalInterface>, bool> LogicHandler<Particle_T>::selectConfiguration(
    Functor &functor, const InteractionTypeOption &interactionType) {
  auto &autoTuner = *_autoTunerRefs[interactionType];

  // Todo: Make LiveInfo persistent between multiple functor calls in the same timestep (e.g. 2B + 3B)
  // https://github.com/AutoPas/AutoPas/issues/916
  LiveInfo info{};
#ifdef AUTOPAS_LOG_LIVEINFO
  auto particleIter = this->begin(IteratorBehavior::ownedOrHalo);
  info.gather(particleIter, _neighborListRebuildFrequency, getNumberOfParticlesOwned(), _logicHandlerInfo.boxMin,
              _logicHandlerInfo.boxMax, _logicHandlerInfo.cutoff, _logicHandlerInfo.verletSkin);
  _liveInfoLogger.logLiveInfo(info, _iteration);
#endif

  // if this iteration is not relevant, take the same algorithm config as before.
  if (not functor.isRelevantForTuning()) {
    auto configuration = autoTuner.getCurrentConfig();
    auto [traversalPtr, _] = isConfigurationApplicable(configuration, functor);

    if (not traversalPtr) {
      // TODO: Can we handle this case gracefully?
      utils::ExceptionHandler::exception(
          "LogicHandler: Functor {} is not relevant for tuning but the given configuration is not applicable!",
          functor.getName());
    }
    return {configuration, std::move(traversalPtr), false};
  }

  if (autoTuner.needsLiveInfo()) {
    // If live info has not been gathered yet, gather it now and send it to the tuner.
    if (info.get().empty()) {
      auto particleIter = this->begin(IteratorBehavior::ownedOrHalo);
      info.gather(particleIter, _neighborListRebuildFrequency, getNumberOfParticlesOwned(), _logicHandlerInfo.boxMin,
                  _logicHandlerInfo.boxMax, _logicHandlerInfo.cutoff, _logicHandlerInfo.verletSkin);
    }
    autoTuner.receiveLiveInfo(info);
  }

  size_t numRejectedConfigs = 0;
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  utils::Timer selectConfigurationTimer;
  selectConfigurationTimer.start();
#endif
  auto [configuration, stillTuning] = autoTuner.getNextConfig();

  // loop as long as we don't get a valid configuration
  do {
    // applicability check also sets the container
    auto [traversalPtr, rejectIndefinitely] = isConfigurationApplicable(configuration, functor);
    if (traversalPtr) {
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
      selectConfigurationTimer.stop();
      AutoPasLog(TRACE, "Select Configuration took {} ms. A total of {} configurations were rejected.",
                 selectConfigurationTimer.getTotalTime(), numRejectedConfigs);
#endif
      return {configuration, std::move(traversalPtr), stillTuning};
    }
    // if no config is left after rejecting this one, an exception is thrown here.
    numRejectedConfigs++;
    std::tie(configuration, stillTuning) = autoTuner.rejectConfig(configuration, rejectIndefinitely);
  } while (true);
}

template <typename Particle_T>
void LogicHandler<Particle_T>::setCurrentContainer(
    std::unique_ptr<ParticleContainerInterface<Particle_T>> newContainer) {
  // copy particles so they do not get lost when the container is switched
  if (_currentContainer != nullptr and newContainer != nullptr) {
    // with these assumptions slightly more space is reserved as numParticlesTotal already includes halos
    const auto numParticlesTotal = _currentContainer->size();
    const auto numParticlesHalo = utils::NumParticlesEstimator::estimateNumHalosUniform(
        numParticlesTotal, _currentContainer->getBoxMin(), _currentContainer->getBoxMax(),
        _currentContainer->getInteractionLength());

    newContainer->reserve(numParticlesTotal, numParticlesHalo);
    for (auto particleIter = _currentContainer->begin(IteratorBehavior::ownedOrHalo); particleIter.isValid();
         ++particleIter) {
      // add a particle as inner if it is owned
      if (particleIter->isOwned()) {
        newContainer->addParticle(*particleIter);
      } else {
        newContainer->addHaloParticle(*particleIter);
      }
    }
  }

  _currentContainer = std::move(newContainer);
}

template <typename Particle_T>
template <class Functor>
bool LogicHandler<Particle_T>::computeInteractionsPipeline(Functor *functor,
                                                           const InteractionTypeOption &interactionType) {
  if (not _interactionTypes.count(interactionType)) {
    utils::ExceptionHandler::exception(
        "LogicHandler::computeInteractionsPipeline(): AutPas was not initialized for the Functor's interactions type: "
        "{}.",
        interactionType);
  }
  /// Selection of configuration (tuning if necessary)
  utils::Timer tuningTimer;
  tuningTimer.start();
  const auto [configuration, traversalPtr, stillTuning] = selectConfiguration(*functor, interactionType);
  tuningTimer.stop();
  auto &autoTuner = *_autoTunerRefs[interactionType];
  autoTuner.logTuningResult(stillTuning, tuningTimer.getTotalTime());

  // Retrieve rebuild info before calling `computeInteractions()` to get the correct value.
  const auto rebuildIteration = not _neighborListsAreValid.load(std::memory_order_relaxed);

  /// Computing the particle interactions
  AutoPasLog(DEBUG, "Iterating with configuration: {} tuning: {}", configuration.toString(), stillTuning);
  const IterationMeasurements measurements = computeInteractions(*functor, *traversalPtr);

  /// Debug Output
  auto bufferSizeListing = [](const auto &buffers) -> std::string {
    std::stringstream ss;
    size_t sum = 0;
    for (const auto &buffer : buffers) {
      ss << buffer.size() << ", ";
      sum += buffer.size();
    }
    ss << " Total: " << sum;
    return ss.str();
  };
  AutoPasLog(TRACE, "particleBuffer     size : {}", bufferSizeListing(_particleBuffer));
  AutoPasLog(TRACE, "haloParticleBuffer size : {}", bufferSizeListing(_haloParticleBuffer));
  AutoPasLog(DEBUG, "Type of interaction :          {}", interactionType.to_string());
  AutoPasLog(DEBUG, "Container::computeInteractions took {} ns", measurements.timeComputeInteractions);
  AutoPasLog(DEBUG, "RemainderTraversal             took {} ns", measurements.timeRemainderTraversal);
  AutoPasLog(DEBUG, "RebuildNeighborLists           took {} ns", measurements.timeRebuild);
  AutoPasLog(DEBUG, "AutoPas::computeInteractions   took {} ns", measurements.timeTotal);
  if (measurements.energyMeasurementsPossible) {
    AutoPasLog(DEBUG, "Energy Consumption: Watts: {} Joules: {} Seconds: {}", measurements.energyWatts,
               measurements.energyJoules, measurements.energyDeltaT);
  }
  _iterationLogger.logIteration(configuration, _iteration, functor->getName(), stillTuning, tuningTimer.getTotalTime(),
                                measurements);

  _flopLogger.logIteration(_iteration, functor->getNumFLOPs(), functor->getHitRate());

  /// Pass on measurements
  // if this was a major iteration add measurements
  if (functor->isRelevantForTuning()) {
    if (stillTuning) {
      // choose the metric of interest
      const auto measurement = [&]() {
        switch (autoTuner.getTuningMetric()) {
          case TuningMetricOption::time:
            return measurements.timeTotal;
          case TuningMetricOption::energy:
            return measurements.energyTotal;
          default:
            utils::ExceptionHandler::exception("LogicHandler::computeInteractionsPipeline(): Unknown tuning metric.");
            return 0l;
        }
      }();
      autoTuner.addMeasurement(measurement, rebuildIteration);
    }
  } else {
    AutoPasLog(TRACE, "Skipping adding of sample because functor is not marked relevant.");
  }
  ++_functorCalls;
  return stillTuning;
}

template <typename Particle_T>
template <class Functor>
std::tuple<std::unique_ptr<TraversalInterface>, bool> LogicHandler<Particle_T>::isConfigurationApplicable(
    const Configuration &config, Functor &functor) {
  // Check if the container supports the traversal
  const auto allContainerTraversals =
      compatibleTraversals::allCompatibleTraversals(config.container, config.interactionType);
  if (allContainerTraversals.find(config.traversal) == allContainerTraversals.end()) {
    AutoPasLog(WARN, "Configuration rejected: Container {} does not support the traversal {}.", config.container,
               config.traversal);
    return {nullptr, /*rejectIndefinitely*/ true};
  }

  // Check if the functor supports the required Newton 3 mode
  if ((config.newton3 == Newton3Option::enabled and not functor.allowsNewton3()) or
      (config.newton3 == Newton3Option::disabled and not functor.allowsNonNewton3())) {
    AutoPasLog(DEBUG, "Configuration rejected: The functor doesn't support Newton 3 {}!", config.newton3);
    return {nullptr, /*rejectIndefinitely*/ true};
  }

  std::unique_ptr<ParticleContainerInterface<Particle_T>> containerPtr{nullptr};
  auto containerInfo =
      ContainerSelectorInfo(_currentContainer->getBoxMin(), _currentContainer->getBoxMax(),
                            _currentContainer->getCutoff(), config.cellSizeFactor, _currentContainer->getVerletSkin(),
                            _verletClusterSize, _sortingThreshold, config.loadEstimator);

  // If we have no current container or needs to be updated to the new config.container, we need to generate a new
  // container.
  const bool generateNewContainer = _currentContainer == nullptr or
                                    _currentContainer->getContainerType() != config.container or
                                    containerInfo != _currentContainerSelectorInfo;

  if (generateNewContainer) {
    // For now, set the local containerPtr to the new container. We do not copy the particles over and set the member
    // _currentContainer until after we know that the traversal is applicable to the domain.
    containerPtr = ContainerSelector<Particle_T>::generateContainer(config.container, containerInfo);
  }

  const auto traversalInfo =
      generateNewContainer ? containerPtr->getTraversalSelectorInfo() : _currentContainer->getTraversalSelectorInfo();

  // Generates a traversal if applicable, otherwise returns a nullptr
  auto traversalPtr =
      TraversalSelector::generateTraversalFromConfig<Particle_T, Functor>(config, functor, traversalInfo);

  // If the traversal is applicable to the domain, and the configuration requires generating a new container,
  // update the member _currentContainer with setCurrentContainer, copying the particle data over, and update
  // _currentContainerSelectorInfo.
  if (traversalPtr and generateNewContainer) {
    _currentContainerSelectorInfo = containerInfo;
    setCurrentContainer(std::move(containerPtr));
  }

  return {std::move(traversalPtr), /*rejectIndefinitely*/ false};
}

}  // namespace autopas
