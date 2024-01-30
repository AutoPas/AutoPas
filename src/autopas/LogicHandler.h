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
#include "autopas/particles/Particle.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/ContainerSelectorInfo.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/NumParticlesEstimator.h"
#include "autopas/utils/SimilarityFunctions.h"
#include "autopas/utils/StaticCellSelector.h"
#include "autopas/utils/StaticContainerSelector.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/logging/IterationLogger.h"
#include "autopas/utils/logging/Logger.h"
#include "autopas/utils/markParticleAsDeleted.h"

namespace autopas {

/**
 * The LogicHandler takes care of the containers s.t. they are all in the same valid state.
 * This is mainly done by incorporating a global container rebuild frequency, which defines when containers and their
 * neighbor lists will be rebuild.
 */
template <typename Particle>
class LogicHandler {
 public:
  /**
   * Constructor of the LogicHandler.
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
        _containerSelector(logicHandlerInfo.boxMin, logicHandlerInfo.boxMax, logicHandlerInfo.cutoff),
        _verletClusterSize(logicHandlerInfo.verletClusterSize),
        _sortingThreshold(logicHandlerInfo.sortingThreshold),
        _iterationLogger(outputSuffix),
        _bufferLocks(std::max(2, autopas::autopas_get_max_threads())) {
    using namespace autopas::utils::ArrayMath::literals;

    // Initialize AutoPas with tuners for given interaction types
    for (auto &[interactionType, tuner] : autotuners) {
      _interactionTypes.insert(interactionType);
      _synchronizer.addInteractionType(interactionType);

      const auto configuration = tuner->getCurrentConfig();
      // initialize the container and make sure it is valid
      const ContainerSelectorInfo containerSelectorInfo{
          configuration.cellSizeFactor, _logicHandlerInfo.verletSkinPerTimestep, _neighborListRebuildFrequency,
          _verletClusterSize, configuration.loadEstimator};
      _containerSelector.selectContainer(configuration.container, containerSelectorInfo);
      checkMinimalSize();
    }

    // initialize locks needed for remainder traversal
    const auto boxLength = logicHandlerInfo.boxMax - logicHandlerInfo.boxMin;
    const auto interactionLengthInv =
        1. / (logicHandlerInfo.cutoff + logicHandlerInfo.verletSkinPerTimestep * rebuildFrequency);
    initSpacialLocks(boxLength, interactionLengthInv);
    for (auto &lockPtr : _bufferLocks) {
      lockPtr = std::make_unique<std::mutex>();
    }
  }

  /**
   * Returns a non-const reference to the currently selected particle container.
   * @return Non-const reference to the container.
   */
  inline autopas::ParticleContainerInterface<Particle> &getContainer() {
    return _containerSelector.getCurrentContainer();
  }

  /**
   * Collects leaving particles from buffer and potentially inserts owned particles to the container.
   * @param insertOwnedParticlesToContainer Decides whether to insert owned particles to the container.
   * @return Leaving particles.
   */
  [[nodiscard]] std::vector<Particle> collectLeavingParticlesFromBuffer(bool insertOwnedParticlesToContainer) {
    const auto &boxMin = _containerSelector.getCurrentContainer().getBoxMin();
    const auto &boxMax = _containerSelector.getCurrentContainer().getBoxMax();
    std::vector<Particle> leavingBufferParticles{};
    for (auto &cell : _particleBuffer) {
      auto &buffer = cell._particles;
      if (insertOwnedParticlesToContainer) {
        for (const auto &p : buffer) {
          if (p.isDummy()) {
            continue;
          }
          if (utils::inBox(p.getR(), boxMin, boxMax)) {
            _containerSelector.getCurrentContainer().addParticle(p);
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
            // Do not increment the iter afterwards!
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
  [[nodiscard]] std::vector<Particle> updateContainer() {
    bool doDataStructureUpdate = not neighborListsAreValid();

    if (doDataStructureUpdate) {
      _neighborListsAreValid.store(false, std::memory_order_relaxed);
    }
    // The next call also adds particles to the container if doDataStructureUpdate is true.
    auto leavingBufferParticles = collectLeavingParticlesFromBuffer(doDataStructureUpdate);

    AutoPasLog(DEBUG, "Initiating container update.");
    auto leavingParticles = _containerSelector.getCurrentContainer().updateContainer(not doDataStructureUpdate);
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
  std::vector<Particle> resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
    using namespace autopas::utils::ArrayMath::literals;
    const auto &oldMin = _containerSelector.getCurrentContainer().getBoxMin();
    const auto &oldMax = _containerSelector.getCurrentContainer().getBoxMax();

    // if nothing changed do nothing
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

    // check all particles
    std::vector<Particle> particlesNowOutside;
    for (auto pIter = _containerSelector.getCurrentContainer().begin(); pIter.isValid(); ++pIter) {
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

    // actually resize the container
    _containerSelector.resizeBox(boxMin, boxMax);
    // The container might have changed sufficiently enough so that we need more or less spacial locks
    const auto boxLength = boxMax - boxMin;
    const auto interactionLengthInv = 1. / (_containerSelector.getCurrentContainer().getInteractionLength());
    initSpacialLocks(boxLength, interactionLengthInv);

    // Set this flag, s.t., the container is rebuilt!
    _neighborListsAreValid.store(false, std::memory_order_relaxed);

    return particlesNowOutside;
  }

  /**
   * Estimates number of halo particles via autopas::utils::NumParticlesEstimator::estimateNumHalosUniform() then
   * calls LogicHandler::reserve(size_t numParticles, size_t numHaloParticles).
   *
   * @param numParticles Total number of owned particles.
   */
  void reserve(size_t numParticles) {
    const auto &container = _containerSelector.getCurrentContainer();
    const auto numParticlesHaloEstimate = autopas::utils::NumParticlesEstimator::estimateNumHalosUniform(
        numParticles, container.getBoxMin(), container.getBoxMax(), container.getInteractionLength());
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
    // there is currently no good heuristic for this buffer so reuse the one for halos.
    for (auto &buffer : _particleBuffer) {
      buffer.reserve(numHaloParticlesPerBuffer);
    }

    _containerSelector.getCurrentContainer().reserve(numParticles, numHaloParticles);
  }

  /**
   * @copydoc AutoPas::addParticle()
   */
  void addParticle(const Particle &p) {
    // first check that the particle actually belongs in the container
    const auto &boxMin = _containerSelector.getCurrentContainer().getBoxMin();
    const auto &boxMax = _containerSelector.getCurrentContainer().getBoxMax();
    if (utils::notInBox(p.getR(), boxMin, boxMax)) {
      autopas::utils::ExceptionHandler::exception(
          "LogicHandler: Trying to add a particle that is not in the bounding box.\n"
          "Box Min {}\n"
          "Box Max {}\n"
          "{}",
          boxMin, boxMax, p.toString());
    }
    if (not neighborListsAreValid()) {
      // Container has to (about to) be invalid to be able to add Particles!
      _containerSelector.getCurrentContainer().template addParticle<false>(p);
    } else {
      // If the container is valid, we add it to the particle buffer.
      _particleBuffer[autopas_get_thread_num()].addParticle(p);
    }
    _numParticlesOwned.fetch_add(1, std::memory_order_relaxed);
  }

  /**
   * @copydoc AutoPas::addHaloParticle()
   */
  void addHaloParticle(const Particle &haloParticle) {
    auto &container = _containerSelector.getCurrentContainer();
    const auto &boxMin = container.getBoxMin();
    const auto &boxMax = container.getBoxMax();
    if (utils::inBox(haloParticle.getR(), boxMin, boxMax)) {
      autopas::utils::ExceptionHandler::exception(
          "LogicHandler: Trying to add a halo particle that is not outside the box of the container.\n"
          "Box Min {}\n"
          "Box Max {}\n"
          "{}",
          utils::ArrayUtils::to_string(boxMin), utils::ArrayUtils::to_string(boxMax), haloParticle.toString());
    }
    if (not neighborListsAreValid()) {
      // If the neighbor lists are not valid, we can add the particle.
      container.template addHaloParticle</* checkInBox */ false>(haloParticle);
    } else {
      // Check if we can update an existing halo(dummy) particle.
      bool updated = container.updateHaloParticle(haloParticle);
      if (not updated) {
        // If we couldn't find an existing particle, add it to the halo particle buffer.
        _haloParticleBuffer[autopas_get_thread_num()].addParticle(haloParticle);
        _haloParticleBuffer[autopas_get_thread_num()]._particles.back().setOwnershipState(OwnershipState::halo);
      }
    }
    _numParticlesHalo.fetch_add(1, std::memory_order_relaxed);
  }

  /**
   * @copydoc AutoPas::deleteAllParticles()
   */
  void deleteAllParticles() {
    _neighborListsAreValid.store(false, std::memory_order_relaxed);
    _containerSelector.getCurrentContainer().deleteAllParticles();
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
  std::tuple<bool, bool> deleteParticleFromBuffers(Particle &particle) {
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
  void decreaseParticleCounter(Particle &particle) {
    if (particle.isOwned()) {
      _numParticlesOwned.fetch_sub(1, std::memory_order_relaxed);
    } else {
      _numParticlesHalo.fetch_sub(1, std::memory_order_relaxed);
    }
  }

  /**
   * This function covers the full pipeline of all mechanics happening during the pairwise iteration.
   * This includes:
   * - selecting a configuration
   *   - gather live info, homogeneity, and max density
   *   - get next config (tuning)
   *   - check applicability
   *   - instantiation of traversal and container
   * - triggering iteration and tuning result logger
   * - pairwise iteration
   *   - init and end traversal
   *   - remainder traversal
   *   - measurements
   * - pass measurements to tuner
   *
   * @tparam Functor
   * @param functor
   * @return True if this was a tuning iteration.
   */
  template <class Functor>
  bool iteratePairwisePipeline(Functor *functor);

  /**
   * This function covers the full pipeline of all mechanics happening during the triwise iteration.
   * This includes:
   * - selecting a configuration
   *   - gather live info, homogeneity, and max density
   *   - get next config (tuning)
   *   - check applicability
   *   - instantiation of traversal and container
   * - triggering iteration and tuning result logger
   * - triwise iteration
   *   - init and end traversal
   *   - remainder traversal
   *   - measurements
   * - pass measurements to tuner
   *
   * @tparam Functor
   * @param functor 3-body functor
   * @return True if this was a tuning iteration.
   */
  template <class Functor>
  bool iterateTriwisePipeline(Functor *functor);

  /**
   * Create the additional vectors vector for a given iterator behavior.
   * @tparam Iterator
   * @param behavior
   * @return Vector of pointers to buffer vectors.
   */
  template <class Iterator>
  typename Iterator::ParticleVecType gatherAdditionalVectors(IteratorBehavior behavior) {
    typename Iterator::ParticleVecType additionalVectors;
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
    return additionalVectors;
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ContainerIterator<Particle, true, false> begin(IteratorBehavior behavior) {
    auto additionalVectors = gatherAdditionalVectors<ContainerIterator<Particle, true, false>>(behavior);
    return _containerSelector.getCurrentContainer().begin(behavior, &additionalVectors);
  }

  /**
   * @copydoc AutoPas::begin()
   */
  autopas::ContainerIterator<Particle, false, false> begin(IteratorBehavior behavior) const {
    auto additionalVectors =
        const_cast<LogicHandler *>(this)->gatherAdditionalVectors<ContainerIterator<Particle, false, false>>(behavior);
    return _containerSelector.getCurrentContainer().begin(behavior, &additionalVectors);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ContainerIterator<Particle, true, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                     const std::array<double, 3> &higherCorner,
                                                                     IteratorBehavior behavior) {
    // sanity check: Most of our stuff depends on `inBox` which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner[d] > higherCorner[d]) {
        autopas::utils::ExceptionHandler::exception(
            "Requesting region Iterator where the upper corner is lower than the lower corner!\n"
            "Lower corner: {}\n"
            "Upper corner: {}",
            lowerCorner, higherCorner);
      }
    }

    auto additionalVectors = gatherAdditionalVectors<ContainerIterator<Particle, true, true>>(behavior);
    return _containerSelector.getCurrentContainer().getRegionIterator(lowerCorner, higherCorner, behavior,
                                                                      &additionalVectors);
  }

  /**
   * @copydoc AutoPas::getRegionIterator()
   */
  autopas::ContainerIterator<Particle, false, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                      const std::array<double, 3> &higherCorner,
                                                                      IteratorBehavior behavior) const {
    // sanity check: Most of our stuff depends on `inBox` which does not handle lowerCorner > higherCorner well.
    for (size_t d = 0; d < 3; ++d) {
      if (lowerCorner[d] > higherCorner[d]) {
        autopas::utils::ExceptionHandler::exception(
            "Requesting region Iterator where the upper corner is lower than the lower corner!\n"
            "Lower corner: {}\n"
            "Upper corner: {}",
            lowerCorner, higherCorner);
      }
    }

    auto additionalVectors =
        const_cast<LogicHandler *>(this)->gatherAdditionalVectors<ContainerIterator<Particle, false, true>>(behavior);
    return std::as_const(_containerSelector)
        .getCurrentContainer()
        .getRegionIterator(lowerCorner, higherCorner, behavior, &additionalVectors);
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
  bool checkTuningStates(InteractionTypeOption::Value interactionType) {
    return _synchronizer.checkTuningState(interactionType);
  }

  /**
   * Update the internal iteration counters.
   */
  void bumpIterationCounters() {
    _stepsSinceLastListRebuild++;
    _iteration++;
  }

  /**
   * Checks if the given configuration can be used with the given functor and the current state of the simulation.
   *
   * @note For the checks we need to switch to the container in the config, hece this function can't be const.
   * Also we need to build the traversal, hence, it is returned.
   *
   * @tparam Functor
   * @param conf
   * @param functor
   * @return tuple<optional<Traversal>, rejectIndefinitely> The optional is empty if the configuration is not applicable
   * The bool rejectIndefinitely indicates if the configuration can be completely removed from the search space because
   * it will never be applicable.
   */
  template <InteractionTypeOption::Value interactionType, class Functor>
  [[nodiscard]] std::tuple<std::optional<std::unique_ptr<TraversalInterface<interactionType>>>, bool>
  isConfigurationApplicable(const Configuration &conf, Functor &functor);

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
  void setParticleBuffers(const std::vector<FullParticleCell<Particle>> &particleBuffers,
                          const std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Getter for the particle buffers.
   *
   * @note Intended for tests only.
   *
   * @return tuple of const references to the internal buffers.
   */
  std::tuple<const std::vector<FullParticleCell<Particle>> &, const std::vector<FullParticleCell<Particle>> &>
  getParticleBuffers() const;

 private:
  /**
   * Initialize or update the spacial locks used during the remainder traversal.
   * If the locks are already initialized but the container size changed, surplus locks will
   * be deleted, new locks are allocated and locks that are still necessary are reused.
   *
   * @param boxLength
   * @param interactionLengthInv
   */
  void initSpacialLocks(const std::array<double, 3> &boxLength, double interactionLengthInv) {
    using namespace autopas::utils::ArrayMath::literals;
    using autopas::utils::ArrayMath::ceil;
    using autopas::utils::ArrayUtils::static_cast_copy_array;

    // one lock per interaction length + one for each halo region
    const auto locksPerDim = static_cast_copy_array<size_t>(ceil(boxLength * interactionLengthInv) + 2.);
    _spacialLocks.resize(locksPerDim[0]);
    for (auto &lockVecVec : _spacialLocks) {
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
   * @return
   */
  template <InteractionTypeOption::Value interactionType, class Functor>
  std::tuple<Configuration, std::unique_ptr<TraversalInterface<interactionType>>, bool> selectConfiguration(
      Functor &functor);

  /**
   * Helper struct collecting all sorts of measurements taken during the pairwise iteration.
   */
  struct IterationMeasurements {
    /// Time
    long timeIteratePairwise{};
    long timeRemainderTraversal{};
    long timeRebuild{};
    long timeTotal{};
    /// Energy. See RaplMeter.h for the meaning of each field.
    bool energyMeasurementsPossible{false};
    double energyPsys{};
    double energyPkg{};
    double energyRam{};
    long energyTotal{};
  };

  /**
   * Triggers the core steps of the pairwise iteration:
   *    - functor init- / end traversal
   *    - rebuilding of neighbor lists
   *    - container.computeInteractions()
   *    - remainder traversal
   *    - time and energy measurements.
   *
   * @tparam PairwiseFunctor
   * @param functor
   * @param traversal
   * @return Struct containing time and energy measurements. If no energy measurements were possible the respective
   * fields are filled with NaN.
   */
  template <class PairwiseFunctor>
  IterationMeasurements iteratePairwise(PairwiseFunctor &functor,
                                        TraversalInterface<InteractionTypeOption::pairwise> &traversal);

  /**
   * Performs the interactions ParticleContainer::iteratePairwise() did not cover.
   *
   * These interactions are:
   *  - particleBuffer    <-> container
   *  - haloParticleBuffer -> container
   *  - particleBuffer    <-> particleBuffer
   *  - haloParticleBuffer -> particleBuffer
   *
   * @note Buffers need to have at least one (empty) cell. They must not be empty.
   *
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam PairwiseFunctor
   * @param f
   * @param container Reference the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class ContainerType, class PairwiseFunctor>
  void doRemainderTraversal(PairwiseFunctor *f, ContainerType &container,
                            std::vector<FullParticleCell<Particle>> &particleBuffers,
                            std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Helper Method for doRemainderTraversal. This method calculates all interactions between buffers and containers
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam PairwiseFunctor
   * @param f
   * @param container Reference the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class ContainerType, class PairwiseFunctor>
  void remainderHelperBufferContainer(PairwiseFunctor *f, ContainerType &container,
                                      std::vector<FullParticleCell<Particle>> &particleBuffers,
                                      std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Helper Method for doRemainderTraversal. This method calculates all interactions between buffers and buffers
   * @tparam newton3
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class PairwiseFunctor>
  void remainderHelperBufferBuffer(PairwiseFunctor *f, std::vector<FullParticleCell<Particle>> &particleBuffers,
                                   std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Helper Method for doRemainderTraversal. This method calculates all interactions between buffers and halo buffers
   * @tparam newton3
   * @tparam PairwiseFunctor
   * @param f
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class PairwiseFunctor>
  void remainderHelperBufferHaloBuffer(PairwiseFunctor *f, std::vector<FullParticleCell<Particle>> &particleBuffers,
                                       std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Triggers the core steps of the triwise iteration:
   *    - functor init- / end traversal
   *    - rebuilding of neighbor lists
   *    - container.iterateTriwise()
   *    - remainder traversal
   *    - time and energy measurements.
   *
   * @tparam TriwiseFunctor
   * @param functor
   * @param traversal
   * @return Struct containing time and energy measurements. If no energy measurements were possible the respective
   * fields are filled with NaN.
   */
  template <class TriwiseFunctor>
  IterationMeasurements iterateTriwise(TriwiseFunctor &functor,
                                       TraversalInterface<InteractionTypeOption::threeBody> &traversal);

  /**
   * Performs the interactions ParticleContainer::iterateTriwise() did not cover.
   *
   * These interactions are:
   *  - particleBuffer    <-> container
   *  - haloParticleBuffer -> container
   *  - particleBuffer    <-> particleBuffer
   *  - haloParticleBuffer -> particleBuffer
   *
   * @note Buffers need to have at least one (empty) cell. They must not be empty.
   * @note TODO: Update this function with SoA-usage
   *
   * @tparam newton3
   * @tparam ContainerType Type of the particle container.
   * @tparam TriwiseFunctor
   * @param f
   * @param container Reference to the container. Preferably pass the container with the actual static type.
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class ContainerType, class TriwiseFunctor>
  void doRemainderTraversal3B(TriwiseFunctor *f, ContainerType &container,
                              std::vector<FullParticleCell<Particle>> &particleBuffers,
                              std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Check that the simulation box is at least of interaction length in each direction.
   *
   * Throws an exception if minimal size requirements are violated.
   */
  void checkMinimalSize() const;

  /**
   * Checks if in the next iteration the neighbor lists have to be rebuilt.
   *
   * This can be the case either because we hit the rebuild frequency or because the auto tuner tests
   * a new configuration.
   *
   * @return True iff the neighbor lists will not be rebuild.
   */
  bool neighborListsAreValid();

  const LogicHandlerInfo &_logicHandlerInfo;
  /**
   * Specifies after how many pair-wise traversals the neighbor lists (if they exist) are to be rebuild.
   */
  unsigned int _neighborListRebuildFrequency;

  /**
   * Number of particles in a VCL cluster.
   */
  unsigned int _verletClusterSize;

  /**
   * Number of particles in two cells from which sorting should be performed for traversal that use the CellFunctor
   */
  size_t _sortingThreshold;

  /**
   * Reference to the map of AutoTuners held by AutoPas.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> &_autoTunerRefs;

  std::set<InteractionTypeOption> _interactionTypes{};

  autopas::TunerSynchronizer _synchronizer{};

  /**
   * Specifies if the neighbor list is valid.
   */
  std::atomic<bool> _neighborListsAreValid{false};

  /**
   * Steps since last rebuild
   */
  unsigned int _stepsSinceLastListRebuild{std::numeric_limits<unsigned int>::max()};

  /**
   * The current iteration number.
   */
  unsigned int _iteration{0};

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
  std::vector<FullParticleCell<Particle>> _particleBuffer;

  /**
   * Buffer to store halo particles that should not yet be added to the container. There is one buffer per thread.
   */
  std::vector<FullParticleCell<Particle>> _haloParticleBuffer;

  /**
   * Object holding the actual particle container with the ability to switch container types.
   */
  ContainerSelector<Particle> _containerSelector;

  /**
   * Locks for regions in the domain. Used for buffer <-> container interaction.
   */
  std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> _spacialLocks;

  /**
   * Locks for the particle buffers.
   */
  std::vector<std::unique_ptr<std::mutex>> _bufferLocks;

  IterationLogger _iterationLogger;
};

template <typename Particle>
void LogicHandler<Particle>::checkMinimalSize() const {
  const auto &container = _containerSelector.getCurrentContainer();
  // check boxSize at least cutoff + skin
  for (unsigned int dim = 0; dim < 3; ++dim) {
    if (container.getBoxMax()[dim] - container.getBoxMin()[dim] < container.getInteractionLength()) {
      autopas::utils::ExceptionHandler::exception(
          "Box (boxMin[{}]={} and boxMax[{}]={}) is too small.\nHas to be at least cutoff({}) + skin({}) = {}.", dim,
          container.getBoxMin()[dim], dim, container.getBoxMax()[dim], container.getCutoff(), container.getVerletSkin(),
          container.getCutoff() + container.getVerletSkin());
    }
  }
}

template <typename Particle>
bool LogicHandler<Particle>::neighborListsAreValid() {
  // TODO: might need to be separated for 3-body - maybe move logic to AutoTuner
  auto needPairRebuild = _interactionTypes.count(InteractionTypeOption::pairwise) != 0 &&
                         _autoTunerRefs[InteractionTypeOption::pairwise]->willRebuildNeighborLists();
  auto needTriRebuild = _interactionTypes.count(InteractionTypeOption::threeBody) != 0 &&
                        _autoTunerRefs[InteractionTypeOption::threeBody]->willRebuildNeighborLists();

  if (_stepsSinceLastListRebuild >= _neighborListRebuildFrequency or needPairRebuild or needTriRebuild) {
    _neighborListsAreValid.store(false, std::memory_order_relaxed);
  }

  return _neighborListsAreValid.load(std::memory_order_relaxed);
}

template <typename Particle>
void LogicHandler<Particle>::setParticleBuffers(const std::vector<FullParticleCell<Particle>> &particleBuffers,
                                                const std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
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

template <typename Particle>
std::tuple<const std::vector<FullParticleCell<Particle>> &, const std::vector<FullParticleCell<Particle>> &>
LogicHandler<Particle>::getParticleBuffers() const {
  return {_particleBuffer, _haloParticleBuffer};
}

template <typename Particle>
template <class PairwiseFunctor>
typename LogicHandler<Particle>::IterationMeasurements LogicHandler<Particle>::iteratePairwise(
    PairwiseFunctor &functor, TraversalInterface<InteractionTypeOption::pairwise> &traversal) {
  auto &autoTuner = _autoTunerRefs[InteractionTypeOption::pairwise];
  const bool doListRebuild = not neighborListsAreValid();
  const auto &configuration = autoTuner->getCurrentConfig();
  auto &container = _containerSelector.getCurrentContainer();

  autopas::utils::Timer timerTotal;
  autopas::utils::Timer timerRebuild;
  autopas::utils::Timer timerIteratePairwise;
  autopas::utils::Timer timerRemainderTraversal;

  const bool energyMeasurementsPossible = autoTuner->resetEnergy();
  timerTotal.start();

  functor.initTraversal();
  if (doListRebuild) {
    timerRebuild.start();
    container.rebuildNeighborLists(&traversal);
    timerRebuild.stop();
  }
  timerIteratePairwise.start();
  container.iteratePairwise(&traversal);
  timerIteratePairwise.stop();

  timerRemainderTraversal.start();
  withStaticContainerType(container, [&](auto &actualContainerType) {
    if (configuration.newton3) {
      doRemainderTraversal<true>(&functor, actualContainerType, _particleBuffer, _haloParticleBuffer);
    } else {
      doRemainderTraversal<false>(&functor, actualContainerType, _particleBuffer, _haloParticleBuffer);
    }
  });
  timerRemainderTraversal.stop();
  functor.endTraversal(configuration.newton3);

  const auto [energyPsys, energyPkg, energyRam, energyTotal] = autoTuner->sampleEnergy();

  timerTotal.stop();

  constexpr auto nanD = std::numeric_limits<double>::quiet_NaN();
  constexpr auto nanL = std::numeric_limits<long>::quiet_NaN();
  return {timerIteratePairwise.getTotalTime(),
          timerRemainderTraversal.getTotalTime(),
          timerRebuild.getTotalTime(),
          timerTotal.getTotalTime(),
          energyMeasurementsPossible,
          energyMeasurementsPossible ? energyPsys : nanD,
          energyMeasurementsPossible ? energyPkg : nanD,
          energyMeasurementsPossible ? energyRam : nanD,
          energyMeasurementsPossible ? energyTotal : nanL};
}

template <class Particle>
template <bool newton3, class ContainerType, class PairwiseFunctor>
void LogicHandler<Particle>::doRemainderTraversal(PairwiseFunctor *f, ContainerType &container,
                                                  std::vector<FullParticleCell<Particle>> &particleBuffers,
                                                  std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  // Sanity check. If this is violated feel free to add some logic here that adapts the number of locks.
  if (_bufferLocks.size() < particleBuffers.size()) {
    utils::ExceptionHandler::exception("Not enough locks for non-halo buffers! Num Locks: {}, Buffers: {}",
                                       _bufferLocks.size(), particleBuffers.size());
  }

  // Balance buffers. This makes processing them with static scheduling quite efficient.
  // Also, if particles were not inserted in parallel, this enables us to process them in parallel now.
  // Cost is at max O(2N) worst O(N) per buffer collection and negligible compared to interacting them.
  auto cellToVec = [](auto &cell) -> std::vector<Particle> & { return cell._particles; };
  utils::ArrayUtils::balanceVectors(particleBuffers, cellToVec);
  utils::ArrayUtils::balanceVectors(haloParticleBuffers, cellToVec);

  // The following part performs the main remainder traversal. The actual calculation is done in 4 steps carried out
  // in three helper functions.

  // only activate time measurements if it will actually be logged
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  autopas::utils::Timer timerBufferContainer;
  autopas::utils::Timer timerPBufferPBuffer;
  autopas::utils::Timer timerPBufferHBuffer;

  timerBufferContainer.start();
#endif
  // steps 1 & 2. particleBuffer with all close particles in container and haloParticleBuffer with owned, close
  // particles in container
  remainderHelperBufferContainer<newton3>(f, container, particleBuffers, haloParticleBuffers);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerBufferContainer.stop();
  timerPBufferPBuffer.start();
#endif

  // step 3. particleBuffer with itself and all other buffers
  remainderHelperBufferBuffer<newton3>(f, particleBuffers, haloParticleBuffers);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerPBufferPBuffer.stop();
  timerPBufferHBuffer.start();
#endif

  // step 4. particleBuffer with haloParticleBuffer
  remainderHelperBufferHaloBuffer<newton3>(f, particleBuffers, haloParticleBuffers);

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerPBufferHBuffer.stop();
#endif

  // unpack particle SoAs. Halo data is not interesting
  for (auto &buffer : particleBuffers) {
    f->SoAExtractor(buffer, buffer._particleSoABuffer, 0);
  }

  AutoPasLog(TRACE, "Timer Buffers <-> Container (1+2): {}", timerBufferContainer.getTotalTime());
  AutoPasLog(TRACE, "Timer PBuffers<-> PBuffer   (  3): {}", timerPBufferPBuffer.getTotalTime());
  AutoPasLog(TRACE, "Timer PBuffers<-> HBuffer   (  4): {}", timerPBufferHBuffer.getTotalTime());

  // Note: haloParticleBuffer with itself is NOT needed, as interactions between halo particles are unneeded!
}

template <class Particle>
template <bool newton3, class ContainerType, class PairwiseFunctor>
void LogicHandler<Particle>::remainderHelperBufferContainer(
    PairwiseFunctor *f, ContainerType &container, std::vector<FullParticleCell<Particle>> &particleBuffers,
    std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  using autopas::utils::ArrayUtils::static_cast_copy_array;
  using namespace autopas::utils::ArrayMath::literals;

  const auto haloBoxMin = container.getBoxMin() - container.getInteractionLength();
  const auto interactionLengthInv = 1. / container.getInteractionLength();

  const double cutoff = container.getCutoff();
  // one halo and particle buffer pair per thread
  AUTOPAS_OPENMP(parallel for schedule(static, 1) shared(f, _spacialLocks, haloBoxMin, interactionLengthInv))
  for (int bufferId = 0; bufferId < particleBuffers.size(); ++bufferId) {
    auto &particleBuffer = particleBuffers[bufferId];
    auto &haloParticleBuffer = haloParticleBuffers[bufferId];
    // 1. particleBuffer with all close particles in container
    for (auto &&p1 : particleBuffer) {
      const auto pos = p1.getR();
      const auto min = pos - cutoff;
      const auto max = pos + cutoff;
      container.forEachInRegion(
          [&](auto &p2) {
            const auto lockCoords = static_cast_copy_array<size_t>((p2.getR() - haloBoxMin) * interactionLengthInv);
            if constexpr (newton3) {
              const std::lock_guard<std::mutex> lock(*_spacialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]]);
              f->AoSFunctor(p1, p2, true);
            } else {
              f->AoSFunctor(p1, p2, false);
              // no need to calculate force enacted on a halo
              if (not p2.isHalo()) {
                const std::lock_guard<std::mutex> lock(*_spacialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]]);
                f->AoSFunctor(p2, p1, false);
              }
            }
          },
          min, max, IteratorBehavior::ownedOrHalo);
    }

    // 2. haloParticleBuffer with owned, close particles in container
    for (auto &&p1halo : haloParticleBuffer) {
      const auto pos = p1halo.getR();
      const auto min = pos - cutoff;
      const auto max = pos + cutoff;
      container.forEachInRegion(
          [&](auto &p2) {
            const auto lockCoords = static_cast_copy_array<size_t>((p2.getR() - haloBoxMin) * interactionLengthInv);
            // No need to apply anything to p1halo
            //   -> AoSFunctor(p1, p2, false) not needed as it neither adds force nor Upot (potential energy)
            //   -> newton3 argument needed for correct globals
            const std::lock_guard<std::mutex> lock(*_spacialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]]);
            f->AoSFunctor(p2, p1halo, newton3);
          },
          min, max, IteratorBehavior::owned);
    }
  }
}

template <class Particle>
template <bool newton3, class PairwiseFunctor>
void LogicHandler<Particle>::remainderHelperBufferBuffer(PairwiseFunctor *f,
                                                         std::vector<FullParticleCell<Particle>> &particleBuffers,
                                                         std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  // All (halo-)buffer interactions shall happen vectorized, hence, load all buffer data into SoAs
  for (auto &buffer : particleBuffers) {
    f->SoALoader(buffer, buffer._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }
  for (auto &buffer : haloParticleBuffers) {
    f->SoALoader(buffer, buffer._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }

  AUTOPAS_OPENMP(parallel) {
    // For buffer interactions where bufferA == bufferB we can always enable newton3. For all interactions between
    // different buffers we turn newton3 always off, which ensures that only one thread at a time is writing to a
    // buffer. This saves expensive locks.
    // we can not use collapse here without locks, otherwise races would occur.
    AUTOPAS_OPENMP(for)
    for (size_t i = 0; i < particleBuffers.size(); ++i) {
      for (size_t jj = 0; jj < particleBuffers.size(); ++jj) {
        auto *particleBufferSoAA = &particleBuffers[i]._particleSoABuffer;
        const auto j = (i + jj) % particleBuffers.size();
        if (i == j) {
          f->SoAFunctorSingle(*particleBufferSoAA, true);
        } else {
          auto *particleBufferSoAB = &particleBuffers[j]._particleSoABuffer;
          f->SoAFunctorPair(*particleBufferSoAA, *particleBufferSoAB, false);
        }
      }
    }
  }
}

template <class Particle>
template <bool newton3, class PairwiseFunctor>
void LogicHandler<Particle>::remainderHelperBufferHaloBuffer(
    PairwiseFunctor *f, std::vector<FullParticleCell<Particle>> &particleBuffers,
    std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  // Here, phase / color based parallelism turned out to be more efficient than tasks
  AUTOPAS_OPENMP(parallel)
  for (int interactionOffset = 0; interactionOffset < haloParticleBuffers.size(); ++interactionOffset) {
    AUTOPAS_OPENMP(for)
    for (size_t i = 0; i < particleBuffers.size(); ++i) {
      auto &particleBufferSoA = particleBuffers[i]._particleSoABuffer;
      auto &haloBufferSoA =
          haloParticleBuffers[(i + interactionOffset) % haloParticleBuffers.size()]._particleSoABuffer;
      f->SoAFunctorPair(particleBufferSoA, haloBufferSoA, false);
    }
  }
}

template <typename Particle>
template <class TriwiseFunctor>
typename LogicHandler<Particle>::IterationMeasurements LogicHandler<Particle>::iterateTriwise(
    TriwiseFunctor &functor, TraversalInterface<InteractionTypeOption::threeBody> &traversal) {
  const bool doListRebuild = not neighborListsAreValid();
  const auto &configuration = _autoTunerRefs[InteractionTypeOption::threeBody]->getCurrentConfig();
  auto &container = _containerSelector.getCurrentContainer();

  autopas::utils::Timer timerTotal;
  autopas::utils::Timer timerRebuild;
  autopas::utils::Timer timerIterateTriwise;
  autopas::utils::Timer timerRemainderTraversal;

  const bool energyMeasurementsPossible = _autoTunerRefs[InteractionTypeOption::threeBody]->resetEnergy();
  timerTotal.start();

  functor.initTraversal();
  if (doListRebuild) {
    timerRebuild.start();
    container.rebuildNeighborLists(&traversal);
    timerRebuild.stop();
  }
  timerIterateTriwise.start();
  container.iterateTriwise(&traversal);
  timerIterateTriwise.stop();

  timerRemainderTraversal.start();
  withStaticContainerType(container, [&](auto &actualContainerType) {
    if (configuration.newton3) {
      doRemainderTraversal3B<true>(&functor, actualContainerType, _particleBuffer, _haloParticleBuffer);
    } else {
      doRemainderTraversal3B<false>(&functor, actualContainerType, _particleBuffer, _haloParticleBuffer);
    }
  });
  timerRemainderTraversal.stop();
  functor.endTraversal(configuration.newton3);

  const auto [energyPsys, energyPkg, energyRam, energyTotal] =
      _autoTunerRefs[InteractionTypeOption::threeBody]->sampleEnergy();

  timerTotal.stop();

  constexpr auto nanD = std::numeric_limits<double>::quiet_NaN();
  constexpr auto nanL = std::numeric_limits<long>::quiet_NaN();
  return {timerIterateTriwise.getTotalTime(),
          timerRemainderTraversal.getTotalTime(),
          timerRebuild.getTotalTime(),
          timerTotal.getTotalTime(),
          energyMeasurementsPossible,
          energyMeasurementsPossible ? energyPsys : nanD,
          energyMeasurementsPossible ? energyPkg : nanD,
          energyMeasurementsPossible ? energyRam : nanD,
          energyMeasurementsPossible ? energyTotal : nanL};
}

template <class Particle>
template <bool newton3, class ContainerType, class TriwiseFunctor>
void LogicHandler<Particle>::doRemainderTraversal3B(TriwiseFunctor *f, ContainerType &container,
                                                    std::vector<FullParticleCell<Particle>> &particleBuffers,
                                                    std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  using autopas::utils::ArrayUtils::static_cast_copy_array;
  using namespace autopas::utils::ArrayMath::literals;

  // Vector to collect pointers to all buffer particles
  std::vector<Particle *> bufferParticles;

  auto cellToVec = [](auto &cell) -> std::vector<Particle> & { return cell._particles; };

  const size_t numOwnedBufferParticles =
      std::transform_reduce(particleBuffers.begin(), particleBuffers.end(), 0, std::plus<>(),
                            [&](auto &vec) { return cellToVec(vec).size(); });
  const size_t numHaloBufferParticles =
      std::transform_reduce(haloParticleBuffers.begin(), haloParticleBuffers.end(), 0, std::plus<>(),
                            [&](auto &vec) { return cellToVec(vec).size(); });
  const size_t numTotal = numOwnedBufferParticles + numHaloBufferParticles;

  // Place pointers to all owned and halo particles from the buffers in one vector
  bufferParticles.reserve(numTotal);
  for (auto &buffer : particleBuffers) {
    for (auto &particle : buffer._particles) {
      bufferParticles.push_back(&particle);
    }
  }
  for (auto &buffer : haloParticleBuffers) {
    for (auto &particle : buffer._particles) {
      bufferParticles.push_back(&particle);
    }
  }

  // The following part performs the main remainder traversal.

  // only activate time measurements if it will actually be logged
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  autopas::utils::Timer timerBufferBufferBuffer;
  autopas::utils::Timer timerBufferBufferContainer;
  autopas::utils::Timer timerBufferContainerContainer;

  timerBufferBufferBuffer.start();
#endif
  // Step 1: 3-body interactions of all particles in the buffers (owned and halo)
  AUTOPAS_OPENMP(parallel for)
  for (auto i = 0; i < numOwnedBufferParticles; ++i) {
    Particle &p1 = *bufferParticles[i];

    for (auto j = 0; j < numTotal; ++j) {
      if (i == j) continue;
      Particle &p2 = *bufferParticles[j];

      for (auto k = j + 1; k < numTotal; ++k) {
        if (k == i) continue;
        Particle &p3 = *bufferParticles[k];

        // no newton3 for now due to race conditions
        f->AoSFunctor(p1, p2, p3, false);
      }
    }
  }

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerBufferBufferBuffer.stop();
  timerBufferBufferContainer.start();
#endif

  // Step 2: 3-body interactions of 2 buffer particles with 1 container particle
  const auto haloBoxMin = container.getBoxMin() - container.getInteractionLength();
  const auto interactionLengthInv = 1. / container.getInteractionLength();

  const double cutoff = container.getCutoff();
  AUTOPAS_OPENMP(parallel for)
  for (auto i = 0; i < numTotal; ++i) {
    Particle &p1 = *bufferParticles[i];
    const auto pos = p1.getR();
    const auto min = pos - cutoff;
    const auto max = pos + cutoff;

    for (auto j = 0; j < numTotal; ++j) {
      if (j == i) continue;
      Particle &p2 = *bufferParticles[j];

      container.forEachInRegion(
          [&](auto &p3) {
            const auto lockCoords = static_cast_copy_array<size_t>((p3.getR() - haloBoxMin) * interactionLengthInv);
            // no newton3 here due to race conditions
            if (i < numOwnedBufferParticles) f->AoSFunctor(p1, p2, p3, false);
            // no need to calculate force enacted on a halo
            if (not p3.isHalo() and i < j) {
              const std::lock_guard<std::mutex> lock(*_spacialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]]);
              f->AoSFunctor(p3, p1, p2, false);
            }
          },
          min, max, IteratorBehavior::ownedOrHalo);
    }
  }

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerBufferBufferContainer.stop();
  timerBufferContainerContainer.start();
#endif

  // Step 3: 3-body interactions of 1 buffer particle and 2 container particles
  // todo: parallelize without race conditions
  // AUTOPAS_OPENMP(parallel for shared(bufferParticles))
  for (auto i = 0; i < numTotal; ++i) {
    Particle &p1 = *bufferParticles[i];
    const auto pos = p1.getR();
    const auto boxmin = pos - cutoff;
    const auto boxmax = pos + cutoff;

    auto p2Iter = container.getRegionIterator(
        boxmin, boxmax, IteratorBehavior::ownedOrHalo | IteratorBehavior::forceSequential, nullptr);
    for (; p2Iter.isValid(); ++p2Iter) {
      Particle &p2 = *p2Iter;

      auto p3Iter = p2Iter;
      ++p3Iter;
      for (; p3Iter.isValid(); ++p3Iter) {
        Particle &p3 = *p3Iter;

        if constexpr (newton3) {
          f->AoSFunctor(p1, p2, p3, true);
        } else {
          if (i < numOwnedBufferParticles) {
            f->AoSFunctor(p1, p2, p3, false);
          }
          // no need to calculate force enacted on a halo
          if (p2.isOwned()) {
            f->AoSFunctor(p2, p1, p3, false);
          }
          if (p3.isOwned()) {
            f->AoSFunctor(p3, p1, p2, false);
          }
        }
      }
    }
  }

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerBufferContainerContainer.stop();
#endif

  AutoPasLog(TRACE, "Timer Buffer <-> Buffer <-> Buffer       : {}", timerBufferBufferBuffer.getTotalTime());
  AutoPasLog(TRACE, "Timer Buffer <-> Buffer <-> Container    : {}", timerBufferBufferContainer.getTotalTime());
  AutoPasLog(TRACE, "Timer Buffer <-> Container <-> Container : {}", timerBufferContainerContainer.getTotalTime());
}

template <typename Particle>
template <InteractionTypeOption::Value interactionType, class Functor>
std::tuple<Configuration, std::unique_ptr<TraversalInterface<interactionType>>, bool>
LogicHandler<Particle>::selectConfiguration(Functor &functor) {
  bool stillTuning = false;
  Configuration configuration{};
  std::optional<std::unique_ptr<TraversalInterface<interactionType>>> traversalPtrOpt{};
  AutoTuner *tuner;

  if constexpr (utils::isPairwiseFunctor<Functor>()) {
    tuner = _autoTunerRefs[InteractionTypeOption::pairwise].get();
  } else if constexpr (utils::isTriwiseFunctor<Functor>()) {
    tuner = _autoTunerRefs[InteractionTypeOption::threeBody].get();
  } else {
    return {configuration, std::move(traversalPtrOpt.value()), stillTuning};
  }

  // if this iteration is not relevant take the same algorithm config as before.
  if (not functor.isRelevantForTuning()) {
    stillTuning = false;
    configuration = tuner->getCurrentConfig();
    if (_containerSelector.getCurrentContainer().getContainerType() != configuration.container) {
      _containerSelector.selectContainer(
          configuration.container,
          ContainerSelectorInfo(
              configuration.cellSizeFactor,
              _containerSelector.getCurrentContainer().getVerletSkin() / _neighborListRebuildFrequency,
              _neighborListRebuildFrequency, _verletClusterSize, configuration.loadEstimator));
    }
    const auto &container = _containerSelector.getCurrentContainer();
    traversalPtrOpt = autopas::utils::withStaticCellType<Particle>(
        container.getParticleCellTypeEnum(), [&](const auto &particleCellDummy) -> decltype(traversalPtrOpt) {
          // Can't make this unique_ptr const otherwise we can't move it later.
          auto traversalPtr =
              TraversalSelector<std::decay_t<decltype(particleCellDummy)>, interactionType>::template generateTraversal<
                  Functor>(configuration.traversal, functor, container.getTraversalSelectorInfo(),
                           configuration.dataLayout, configuration.newton3);

          // set sortingThreshold of the traversal if it can be casted to a CellPairTraversal and uses the CellFunctor
          if (auto *cellTraversalPtr =
                  dynamic_cast<autopas::CellTraversal<std::decay_t<decltype(particleCellDummy)>> *>(
                      traversalPtr.get())) {
            cellTraversalPtr->setSortingThreshold(_sortingThreshold);
          }
          if (traversalPtr->isApplicable()) {
            return std::optional{std::move(traversalPtr)};
          } else {
            return std::nullopt;
          }
        });
  } else {
    if (tuner->needsHomogeneityAndMaxDensityBeforePrepare()) {
      utils::Timer timerCalculateHomogeneity;
      timerCalculateHomogeneity.start();
      const auto &container = _containerSelector.getCurrentContainer();
      const auto [homogeneity, maxDensity] =
          autopas::utils::calculateHomogeneityAndMaxDensity(container, container.getBoxMin(), container.getBoxMax());
      timerCalculateHomogeneity.stop();
      tuner->addHomogeneityAndMaxDensity(homogeneity, maxDensity, timerCalculateHomogeneity.getTotalTime());
    }

    const auto needsLiveInfo = tuner->prepareIteration();

    if (needsLiveInfo) {
      LiveInfo info{};
      info.gather(_containerSelector.getCurrentContainer(), functor, _neighborListRebuildFrequency);
      tuner->receiveLiveInfo(info);
    }

    std::tie(configuration, stillTuning) = tuner->getNextConfig();

    // loop as long as we don't get a valid configuration
    bool rejectIndefinitely = false;
    while (true) {
      // applicability check also sets the container
      std::tie(traversalPtrOpt, rejectIndefinitely) =
          isConfigurationApplicable<interactionType>(configuration, functor);
      if (traversalPtrOpt.has_value()) {
        break;
      }
      // if no config is left after rejecting this one an exception is thrown here.
      std::tie(configuration, stillTuning) = tuner->rejectConfig(configuration, rejectIndefinitely);
    }
  }

  // log tuning status for current tuner
  _synchronizer.recordTuningState(interactionType, stillTuning);

  return {configuration, std::move(traversalPtrOpt.value()), stillTuning};
}

template <typename Particle>
template <class Functor>
bool LogicHandler<Particle>::iteratePairwisePipeline(Functor *functor) {
  if (not _interactionTypes.count(InteractionTypeOption::pairwise)) {
    autopas::utils::ExceptionHandler::exception(
        "LogicHandler::iteratePairwisePipeline(): LogicHandler was not initialized for pairwise interactions.");
  }
  /// Selection of configuration (tuning if necessary)
  utils::Timer tuningTimer;
  tuningTimer.start();
  const auto [configuration, traversalPtr, stillTuning] =
      selectConfiguration<InteractionTypeOption::pairwise>(*functor);
  tuningTimer.stop();
  auto &autoTuner = _autoTunerRefs[InteractionTypeOption::pairwise];
  autoTuner->logIteration(configuration, stillTuning, tuningTimer.getTotalTime());

  /// Pairwise iteration
  AutoPasLog(DEBUG, "Iterating with configuration: {} tuning: {}", configuration.toString(), stillTuning);
  const IterationMeasurements measurements = iteratePairwise(*functor, *traversalPtr.get());

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
  AutoPasLog(DEBUG, "Container::iteratePairwise took {} ns", measurements.timeIteratePairwise);
  AutoPasLog(DEBUG, "RemainderTraversal         took {} ns", measurements.timeRemainderTraversal);
  AutoPasLog(DEBUG, "RebuildNeighborLists       took {} ns", measurements.timeRebuild);
  AutoPasLog(DEBUG, "AutoPas::iteratePairwise took {} ns", measurements.timeTotal);
  if (measurements.energyMeasurementsPossible) {
    AutoPasLog(DEBUG, "Energy Consumption: Psys: {} Joules Pkg: {} Joules Ram: {} Joules", measurements.energyPsys,
               measurements.energyPkg, measurements.energyRam);
  }
  _iterationLogger.logIteration(configuration, _iteration, functor->getName(), stillTuning,
                                measurements.timeIteratePairwise, measurements.timeRemainderTraversal,
                                measurements.timeRebuild, measurements.timeTotal, tuningTimer.getTotalTime(),
                                measurements.energyPsys, measurements.energyPkg, measurements.energyRam);

  /// Pass on measurements
  // if this was a major iteration add measurements and bump counters
  if (functor->isRelevantForTuning()) {
    if (stillTuning) {
      switch (autoTuner->getTuningMetric()) {
        case TuningMetricOption::time:
          autoTuner->addMeasurement(measurements.timeTotal, not neighborListsAreValid());
          break;
        case TuningMetricOption::energy:
          autoTuner->addMeasurement(measurements.energyTotal, not neighborListsAreValid());
          break;
      }
    } else {
      AutoPasLog(TRACE, "Skipping adding of sample because functor is not marked relevant.");
    }

    // this function depends on LogicHandler's and the AutoTuner's iteration counters,
    // that should not have been updated yet.
    if (not neighborListsAreValid() /*we have done a rebuild now*/) {
      // list is now valid
      _neighborListsAreValid.store(true, std::memory_order_relaxed);
      _stepsSinceLastListRebuild = 0;
    }
    ++_stepsSinceLastListRebuild;
  }
  return stillTuning;
}

template <typename Particle>
template <class Functor>
bool LogicHandler<Particle>::iterateTriwisePipeline(Functor *functor) {
  if (not _interactionTypes.count(InteractionTypeOption::threeBody)) {
    autopas::utils::ExceptionHandler::exception(
        "LogicHandler::iterateTriwisePipeline(): LogicHandler was not initialized for 3-body interactions.");
  }
  /// Selection of configuration (tuning if necessary)
  utils::Timer tuningTimer;
  tuningTimer.start();
  //  bool stillTuning = true;
  const auto [configuration, traversalPtr, stillTuning] =
      selectConfiguration<InteractionTypeOption::threeBody>(*functor);
  tuningTimer.stop();
  AutoPasLog(DEBUG, "Selecting a configuration took {} ns.", tuningTimer.getTotalTime());
  auto &autoTuner = _autoTunerRefs[InteractionTypeOption::threeBody];
  autoTuner->logIteration(configuration, stillTuning, tuningTimer.getTotalTime());

  /// Triwise iteration
  AutoPasLog(DEBUG, "Iterating with configuration: {} tuning: {}", configuration.toString(), stillTuning);
  const IterationMeasurements measurements = iterateTriwise(*functor, *traversalPtr.get());

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
  AutoPasLog(DEBUG, "Container::iterateTriwise took {} ns", measurements.timeIteratePairwise);
  AutoPasLog(DEBUG, "RemainderTraversal         took {} ns", measurements.timeRemainderTraversal);
  AutoPasLog(DEBUG, "RebuildNeighborLists       took {} ns", measurements.timeRebuild);
  AutoPasLog(DEBUG, "AutoPas::iterateTriwise took {} ns", measurements.timeTotal);
  if (measurements.energyMeasurementsPossible) {
    AutoPasLog(DEBUG, "Energy Consumption: Psys: {} Joules Pkg: {} Joules Ram: {} Joules", measurements.energyPsys,
               measurements.energyPkg, measurements.energyRam);
  }
  _iterationLogger.logIteration(configuration, _iteration, functor->getName(), stillTuning,
                                measurements.timeIteratePairwise, measurements.timeRemainderTraversal,
                                measurements.timeRebuild, measurements.timeTotal, tuningTimer.getTotalTime(),
                                measurements.energyPsys, measurements.energyPkg, measurements.energyRam);

  /// Pass on measurements
  // if this was a major iteration add measurements and bump counters
  if (functor->isRelevantForTuning()) {
    if (stillTuning) {
      switch (autoTuner->getTuningMetric()) {
        case TuningMetricOption::time:
          autoTuner->addMeasurement(measurements.timeTotal, not neighborListsAreValid());
          break;
        case TuningMetricOption::energy:
          autoTuner->addMeasurement(measurements.energyTotal, not neighborListsAreValid());
          break;
      }
    } else {
      AutoPasLog(TRACE, "Skipping adding of sample because functor is not marked relevant.");
    }

    // this function depends on LogicHandler's and the AutoTuner's iteration counters,
    // that should not have been updated yet.
    if (not neighborListsAreValid() /*we have done a rebuild now*/) {
      // list is now valid
      _neighborListsAreValid.store(true, std::memory_order_relaxed);
      _stepsSinceLastListRebuild = 0;
    }
    ++_stepsSinceLastListRebuild;
  }
  return stillTuning;
}

template <typename Particle>
template <InteractionTypeOption::Value interactionType, class Functor>
std::tuple<std::optional<std::unique_ptr<TraversalInterface<interactionType>>>, bool>
LogicHandler<Particle>::isConfigurationApplicable(const Configuration &conf, Functor &functor) {
  // Check if the container supports the traversal
  const auto allContainerTraversals =
      compatibleTraversals::allCompatibleTraversals(conf.container, conf.interactionType);
  if (allContainerTraversals.find(conf.traversal) == allContainerTraversals.end()) {
    AutoPasLog(DEBUG, "Configuration rejected: Container {} does not support the traversal {}.", conf.container,
               conf.traversal);
    return {std::nullopt, true};
  }

  // Check if the required Newton 3 mode is supported by the functor
  if ((conf.newton3 == Newton3Option::enabled and not functor.allowsNewton3()) or
      (conf.newton3 == Newton3Option::disabled and not functor.allowsNonNewton3())) {
    AutoPasLog(DEBUG, "Configuration rejected: The functor doesn't support Newton 3 {}!", conf.newton3);
    return {std::nullopt, true};
  }

  // Check if the traversal is applicable to the current state of the container
  _containerSelector.selectContainer(
      conf.container,
      ContainerSelectorInfo(conf.cellSizeFactor,
                            _containerSelector.getCurrentContainer().getVerletSkin() / _neighborListRebuildFrequency,
                            _neighborListRebuildFrequency, _verletClusterSize, conf.loadEstimator));
  const auto &container = _containerSelector.getCurrentContainer();
  const auto traversalInfo = container.getTraversalSelectorInfo();

  auto traversalPtrOpt = autopas::utils::withStaticCellType<Particle>(
      container.getParticleCellTypeEnum(),
      [&](const auto &particleCellDummy) -> std::optional<std::unique_ptr<TraversalInterface<interactionType>>> {
        // Can't make this unique_ptr const otherwise we can't move it later.
        auto traversalPtr =
            TraversalSelector<std::decay_t<decltype(particleCellDummy)>,
                              interactionType>::template generateTraversal<Functor>(conf.traversal, functor,
                                                                                    traversalInfo, conf.dataLayout,
                                                                                    conf.newton3);

        // set sortingThreshold of the traversal if it can be casted to a CellPairTraversal and uses the CellFunctor
        if (auto *cellTraversalPtr =
                dynamic_cast<autopas::CellTraversal<std::decay_t<decltype(particleCellDummy)>> *>(traversalPtr.get())) {
          cellTraversalPtr->setSortingThreshold(_sortingThreshold);
        }
        if (traversalPtr->isApplicable()) {
          return std::optional{std::move(traversalPtr)};
        } else {
          return std::nullopt;
        }
      });
  return {std::move(traversalPtrOpt), false};
}

}  // namespace autopas
