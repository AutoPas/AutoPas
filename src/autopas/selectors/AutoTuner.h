/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.2018
 */

#pragma once

#include <array>
#include <memory>
#include <set>

#include "Smoothing.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/selectors/TraversalSelector.h"
#include "autopas/selectors/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/StaticCellSelector.h"
#include "autopas/utils/StaticContainerSelector.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/IterationLogger.h"
#include "autopas/utils/logging/TuningDataLogger.h"
#include "autopas/utils/logging/TuningResultLogger.h"

namespace autopas {

/**
 * Calls to the iteratePairwise() method are passed through this class for two reasons:
 * 1. Measuring time of the iteration.
 * 2. Selecting an appropriate configuration for the pairwise iteration.
 *
 * The tuner can be in one of two states. If it currently should look for a new optimum, it is in the
 * so-called tuning phase. During a tuning phase, for each Configuration, multiple measurements can be taken,
 * which are called samples. To reduce noise, the samples for one configuration are then condensed to one value for
 * the current tuning phase, called evidence. The evidences are handed on to a tuningStrategy, which selects a) what
 * Configuration to test next and b) which configuration is the best in this tuning phase.
 * If it should not look for a new optimum it is not in a tuning phase.
 *
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle>
class AutoTuner {
 public:
  /**
   * Constructor for the AutoTuner that generates all configurations from the given options.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff Cutoff radius to be used in this container.
   * @param verletSkinPerTimestep Length added to the cutoff for the Verlet lists' skin.
   * @param verletClusterSize Number of particles in a cluster to use in verlet list.
   * @param tuningStrategy Object implementing the modelling and exploration of a search space.
   * @param MPITuningMaxDifferenceForBucket For MPI-tuning: Maximum of the relative difference in the comparison metric
   * for two ranks which exchange their tuning information.
   * @param MPITuningWeightForMaxDensity For MPI-tuning: Weight for maxDensity in the calculation for bucket
   * distribution.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningInterval Number of time steps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   * @param rebuildFrequency The rebuild frequency this AutoPas instance uses.
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkinPerTimestep,
            unsigned int verletClusterSize, std::unique_ptr<TuningStrategyInterface> tuningStrategy,
            double MPITuningMaxDifferenceForBucket, double MPITuningWeightForMaxDensity,
            SelectorStrategyOption selectorStrategy, unsigned int tuningInterval, unsigned int maxSamples,
            unsigned int rebuildFrequency, const std::string &outputSuffix = "")
      : _selectorStrategy(selectorStrategy),
        _tuningStrategy(std::move(tuningStrategy)),
        _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff),
        _verletSkinPerTimestep(verletSkinPerTimestep),
        _verletClusterSize(verletClusterSize),
        _mpiTuningMaxDifferenceForBucket(MPITuningMaxDifferenceForBucket),
        _mpiTuningWeightForMaxDensity(MPITuningWeightForMaxDensity),
        _rebuildFrequency(rebuildFrequency),
        _maxSamples(maxSamples),
        _samplesNotRebuildingNeighborLists(maxSamples),
        _iteration(0),
        _iterationLogger(outputSuffix),
        _tuningResultLogger(outputSuffix),
        _tuningDataLogger(maxSamples, outputSuffix) {
    using autopas::utils::ArrayMath::ceil;
    using autopas::utils::ArrayMath::mulScalar;
    using autopas::utils::ArrayMath::sub;
    using autopas::utils::ArrayUtils::static_cast_array;

    // initialize locks needed for remainder traversal
    const auto boxLength = sub(boxMax, boxMin);
    const auto interactionLengthInv = 1. / (cutoff + verletSkinPerTimestep * rebuildFrequency);
    const auto locksPerDim = static_cast_array<size_t>(ceil(mulScalar(boxLength, interactionLengthInv)));
    _spacialLocks.resize(locksPerDim[0]);
    for (auto &lockVecVec : _spacialLocks) {
      lockVecVec.resize(locksPerDim[1]);
      for (auto &lockVec : lockVecVec) {
        lockVec.resize(locksPerDim[2]);
        for (auto &lockPtr : lockVec) {
          lockPtr = std::make_unique<std::mutex>();
        }
      }
    }

    if (_tuningStrategy->searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
    }

    selectCurrentContainer();
  }

  /**
   * Move assignment operator
   * @param other
   * @return
   */
  AutoTuner &operator=(AutoTuner &&other) noexcept {
    _tuningStrategy = std::move(other._tuningStrategy);
    return *this;
  }

  /**
   * Pass values to the towards the actual container.
   * @param boxMin
   * @param boxMax
   */
  void resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
    _containerSelector.resizeBox(boxMin, boxMax);
  }

  /**
   * @copydoc autopas::AutoPas::forceRetune()
   */
  void forceRetune();

  /**
   * Getter for the current container.
   * @return Smart pointer to the current container.
   */
  std::shared_ptr<autopas::ParticleContainerInterface<Particle>> getContainer() {
    return _containerSelector.getCurrentContainer();
  }

  /**
   * Getter for the current container.
   * @return Smart pointer to the current container.
   */
  std::shared_ptr<const autopas::ParticleContainerInterface<Particle>> getContainer() const {
    return _containerSelector.getCurrentContainer();
  }

  /**
   * Check if a configuration is applicable to the current domain with the given functor.
   * @tparam PairwiseFunctor
   * @param conf
   * @param pairwiseFunctor
   * @return
   */
  template <class PairwiseFunctor>
  bool configApplicable(const Configuration &conf, PairwiseFunctor &pairwiseFunctor);

  /**
   * Whole tuning logic for one iteration.
   * This function covers:
   *  - Instantiation of the traversal to be used.
   *  - Actual pairwise iteration for application of the functor.
   *  - Management of tuning phases and calls to tune() if necessary.
   *  - Measurement of timings and calls to addTimeMeasurement() if necessary.
   *  - Calls to _iterationLogger if necessary.
   *
   * @note Buffers need to have at least one (empty) cell. They must not be empty.
   *
   * @tparam PairwiseFunctor
   * @param f Functor that describes the pair-potential.
   * @param doListRebuild Indicates whether or not the verlet lists should be rebuild.
   * @param particleBuffers A buffer of additional particles to consider.
   * @param haloParticleBuffers A buffer of additional halo particles to consider.
   * @return true if this was a tuning iteration.
   */
  template <class PairwiseFunctor>
  bool iteratePairwise(PairwiseFunctor *f, bool doListRebuild, std::vector<FullParticleCell<Particle>> &particleBuffers,
                       std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Returns whether the configuration will be changed in the next iteration.
   * This does does not necessarily mean that the container will change.
   *
   * @return True if the next iteratePairwise() call uses a different configuration. False otherwise.
   */
  bool willRebuild() {
    if (_tuningStrategy->searchSpaceIsTrivial()) {
      return false;
    }
    return _iterationsSinceTuning >= _tuningInterval and getCurrentNumSamples() >= _maxSamples;
  }

  /**
   * Get the currently selected configuration.
   * @return
   */
  [[nodiscard]] const Configuration &getCurrentConfig() const;

  /**
   * Initialize the container specified by the TuningStrategy.
   */
  bool searchSpaceIsTrivial();

 private:
  /**
   * Total number of collected samples. This is the sum of the sizes of all sample vectors.
   * @return Sum of sizes of sample vectors.
   */
  auto getCurrentNumSamples() const {
    return _samplesNotRebuildingNeighborLists.size() + _samplesRebuildingNeighborLists.size();
  }

  /**
   * Save the runtime of a given traversal.
   *
   * Samples are collected and reduced to one single value according to _selectorStrategy. Only then the value is passed
   * on to the tuning strategy. This function expects that samples of the same configuration are taken consecutively.
   * The time argument is a long because std::chrono::duration::count returns a long.
   *
   * @param time
   * @param neighborListRebuilt If the neighbor list as been rebuilt during the given time.
   */
  void addTimeMeasurement(long time, bool neighborListRebuilt);

  /**
   * Estimate the runtime from the current samples according to the SelectorStrategy and rebuild frequency.
   * Samples are weighted so that we normalize to the expected number of (non-)rebuild iterations and then divide by the
   * rebuild frequency.
   * @return estimate time for one iteration
   */
  [[nodiscard]] long estimateRuntimeFromSamples() const {
    // reduce samples for rebuild and non-rebuild iterations with the given selector strategy
    const auto reducedValueBuilding =
        autopas::OptimumSelector::optimumValue(_samplesRebuildingNeighborLists, _selectorStrategy);
    // if there is no data for the non rebuild iterations we have to assume them taking the same time as rebuilding ones
    // this might neither be a good estimate nor fair but the best we can do
    const auto reducedValueNotBuilding =
        _samplesNotRebuildingNeighborLists.empty()
            ? reducedValueBuilding
            : autopas::OptimumSelector::optimumValue(_samplesNotRebuildingNeighborLists, _selectorStrategy);

    const auto numIterationsNotBuilding =
        std::max(0, static_cast<int>(_rebuildFrequency) - static_cast<int>(_samplesRebuildingNeighborLists.size()));
    const auto numIterationsBuilding = _rebuildFrequency - numIterationsNotBuilding;

    // calculate weighted estimate for one iteration
    return (numIterationsBuilding * reducedValueBuilding + numIterationsNotBuilding * reducedValueNotBuilding) /
           _rebuildFrequency;
  }

  /**
   * Initialize the container specified by the TuningStrategy.
   */
  void selectCurrentContainer();

  template <class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3, bool inTuningPhase>
  void iteratePairwiseTemplateHelper(PairwiseFunctor *f, bool doListRebuild,
                                     std::vector<FullParticleCell<Particle>> &particleBuffer,
                                     std::vector<FullParticleCell<Particle>> &haloParticleBuffer);

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
   * @tparam T (Smart) pointer Type of the particle container.
   * @tparam PairwiseFunctor
   * @param f
   * @param containerPtr (Smart) Pointer to the container
   * @param particleBuffers vector of particle buffers. These particles' force vectors will be updated.
   * @param haloParticleBuffers vector of halo particle buffers. These particles' force vectors will not necessarily be
   * updated.
   */
  template <bool newton3, class T, class PairwiseFunctor>
  void doRemainderTraversal(PairwiseFunctor *f, T containerPtr,
                            std::vector<FullParticleCell<Particle>> &particleBuffers,
                            std::vector<FullParticleCell<Particle>> &haloParticleBuffers);

  /**
   * Tune available algorithm configurations.
   *
   * It is assumed this function is only called for relevant functors and that at least two configurations are allowed.
   * When in tuning phase selects next config to test. At the end of the tuning phase select optimum.
   * The function returns true if the selected config is not yet the optimum but something that should be sampled.
   *
   * @tparam PairwiseFunctor
   * @param pairwiseFunctor
   * @return true iff still in tuning phase.
   */
  template <class PairwiseFunctor>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  SelectorStrategyOption _selectorStrategy;
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  /**
   * Counter for the simulation iteration.
   */
  size_t _iteration;

  /**
   * Counter for the number of times a tuning phase was started.
   */
  size_t _tuningInterval;

  /**
   * Number of iterations since the end of the last tuning phase.
   */
  size_t _iterationsSinceTuning;

  /**
   * Object holding the actual particle container and having the ability to change it.
   */
  ContainerSelector<Particle> _containerSelector;

  double _verletSkinPerTimestep;
  unsigned int _verletClusterSize;

  /**
   * The rebuild frequency this instance of AutoPas uses.
   */
  unsigned int _rebuildFrequency;

  /**
   * How many times each configuration should be tested.
   */
  const size_t _maxSamples;

  /**
   * Parameter used For MPI-tuning: Maximum of the relative difference in the comparison metric for two ranks which
   * exchange their tuning information.
   */
  double _mpiTuningMaxDifferenceForBucket;

  /**
   * Parameter used For MPI-tuning: Weight for maxDensity in the calculation for bucket distribution.
   */
  double _mpiTuningWeightForMaxDensity;

  /**
   * Buffer for the homogeneities of the last ten Iterations
   */
  std::vector<double> _homogeneitiesOfLastTenIterations;

  /**
   * Buffer for the homogeneities of the last ten Iterations
   */
  std::vector<double> _maxDensitiesOfLastTenIterations;

  /**
   * Raw time samples of the current configuration. Contains only the samples of iterations where the neighbor lists are
   * not rebuilt.
   *
   * @note Initialized with size of _maxSamples to start tuning at start of simulation.
   */
  std::vector<long> _samplesNotRebuildingNeighborLists;

  /**
   * Raw time samples of the current configuration. Contains only the samples of iterations where the neighbor lists
   * have been rebuilt.
   */
  std::vector<long> _samplesRebuildingNeighborLists{};

  /**
   * For each configuration the collection of all evidence (smoothed values) collected so far and in which iteration.
   * Configuration -> vector< iteration, time >
   */
  std::map<Configuration, std::vector<std::pair<size_t, long>>> _evidence;

  /**
   * Timer used to determine how much time is wasted by calculating multiple homogeneities for smoothing
   * Only used temporarily
   */
  autopas::utils::Timer _timerCalculateHomogeneity;

  /**
   * Locks for regions in the domain. Used for buffer <-> container interaction.
   */
  std::vector<std::vector<std::vector<std::unique_ptr<std::mutex>>>> _spacialLocks;

  IterationLogger _iterationLogger;
  TuningResultLogger _tuningResultLogger;
  TuningDataLogger _tuningDataLogger;
};

template <class Particle>
void AutoTuner<Particle>::selectCurrentContainer() {
  auto conf = _tuningStrategy->getCurrentConfiguration();
  _containerSelector.selectContainer(
      conf.container, ContainerSelectorInfo(conf.cellSizeFactor, _verletSkinPerTimestep, _rebuildFrequency,
                                            _verletClusterSize, conf.loadEstimator));
}

/**
 * Access to the searchSpaceIsTrivial bool variable (true if search space size  is 1 or less).
 * @return Smart pointer to the searchSpaceIsTrivial variable.
 */
template <class Particle>
bool AutoTuner<Particle>::searchSpaceIsTrivial() {
  return _tuningStrategy->searchSpaceIsTrivial();
}

template <class Particle>
void AutoTuner<Particle>::forceRetune() {
  _iterationsSinceTuning = _tuningInterval;
  _samplesNotRebuildingNeighborLists.resize(_maxSamples);
}

template <class Particle>
template <class PairwiseFunctor>
bool AutoTuner<Particle>::iteratePairwise(PairwiseFunctor *f, bool doListRebuild,
                                          std::vector<FullParticleCell<Particle>> &particleBuffers,
                                          std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  bool isTuning = false;
  // tune if :
  // - more than one config exists
  // - currently in tuning phase
  // - functor is relevant
  if (_tuningStrategy->smoothedHomogeneityAndMaxDensityNeeded() and _iterationsSinceTuning >= _tuningInterval - 9 and
      _iterationsSinceTuning <= _tuningInterval) {
    _timerCalculateHomogeneity.start();
    const auto &container = getContainer();
    const auto [homogeneity, maxDensity] =
        autopas::utils::calculateHomogeneityAndMaxDensity(*container, container->getBoxMin(), container->getBoxMax());
    _homogeneitiesOfLastTenIterations.push_back(homogeneity);
    _maxDensitiesOfLastTenIterations.push_back(maxDensity);
    _timerCalculateHomogeneity.stop();
  }
  if ((not _tuningStrategy->searchSpaceIsTrivial()) and _iterationsSinceTuning >= _tuningInterval and
      f->isRelevantForTuning()) {
    isTuning = tune<PairwiseFunctor>(*f);
    if (not isTuning) {
      _iterationsSinceTuning = 0;
    }
  }

  // large case differentiation for data layout and newton 3
  switch (_tuningStrategy->getCurrentConfiguration().dataLayout) {
    case DataLayoutOption::aos: {
      if (_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ true,
                                        /*tuning*/ true>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ true,
                                        /*tuning*/ false>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ false,
                                        /*tuning*/ true>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ false,
                                        /*tuning*/ false>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        }
      }
      break;
    }
    case DataLayoutOption::soa: {
      if (_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ true,
                                        /*tuning*/ true>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ true,
                                        /*tuning*/ false>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ false,
                                        /*tuning*/ true>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ false,
                                        /*tuning*/ false>(f, doListRebuild, particleBuffers, haloParticleBuffers);
        }
      }
      break;
    }
    default:
      utils::ExceptionHandler::exception("AutoTuner: Unknown data layout : {}",
                                         _tuningStrategy->getCurrentConfiguration().dataLayout);
  }

  if (f->isRelevantForTuning()) {
    ++_iterationsSinceTuning;
    ++_iteration;
  }
  return isTuning;
}

template <class Particle>
template <bool newton3, class T, class PairwiseFunctor>
void AutoTuner<Particle>::doRemainderTraversal(PairwiseFunctor *f, T containerPtr,
                                               std::vector<FullParticleCell<Particle>> &particleBuffers,
                                               std::vector<FullParticleCell<Particle>> &haloParticleBuffers) {
  using autopas::utils::ArrayMath::addScalar;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::subScalar;
  using autopas::utils::ArrayUtils::static_cast_array;

  const auto boxMin = containerPtr->getBoxMin();
  const auto interactionLengthInv = 1. / containerPtr->getInteractionLength();

  // Balance buffers. This makes processing them with static scheduling quite efficient.
  // Also, if particles were not inserted in parallel, this enables us to process them in parallel now.
  // Cost is at max O(2N) worst O(N) per buffer collection and negligible compared to interacting them.
  auto cellToVec = [](auto &cell) -> std::vector<Particle> & { return cell._particles; };
  utils::ArrayUtils::balanceVectors(particleBuffers, cellToVec);
  utils::ArrayUtils::balanceVectors(haloParticleBuffers, cellToVec);

  // only activate time measurements if it will actually be logged
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  autopas::utils::Timer timerBufferContainer;
  autopas::utils::Timer timerPBufferPBuffer;
  autopas::utils::Timer timerPBufferHBuffer;

  timerBufferContainer.start();
#endif
  withStaticContainerType(containerPtr, [&](auto staticTypedContainerPtr) {
    const double cutoff = staticTypedContainerPtr->getCutoff();
#ifdef AUTOPAS_OPENMP
// one halo and particle buffer pair per thread
#pragma omp parallel for schedule(static, 1), shared(f, _spacialLocks, boxMin, interactionLengthInv)
#endif
    for (int bufferId = 0; bufferId < particleBuffers.size(); ++bufferId) {
      auto &particleBuffer = particleBuffers[bufferId];
      auto &haloParticleBuffer = haloParticleBuffers[bufferId];
      // 1. particleBuffer with all close particles in container
      for (auto &&p1 : particleBuffer) {
        const auto pos = p1.getR();
        const auto min = subScalar(pos, cutoff);
        const auto max = addScalar(pos, cutoff);
        staticTypedContainerPtr->forEachInRegion(
            [&](auto &p2) {
              const auto lockCoords =
                  static_cast_array<size_t>(mulScalar(sub(p2.getR(), boxMin), interactionLengthInv));
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
        const auto min = autopas::utils::ArrayMath::subScalar(pos, cutoff);
        const auto max = autopas::utils::ArrayMath::addScalar(pos, cutoff);
        staticTypedContainerPtr->forEachInRegion(
            [&](auto &p2) {
              const auto lockCoords =
                  static_cast_array<size_t>(mulScalar(sub(p2.getR(), boxMin), interactionLengthInv));
              // No need to apply anything to p1halo
              //   -> AoSFunctor(p1, p2, false) not needed as it neither adds force nor Upot
              //   -> newton3 argument needed for correct globals
              const std::lock_guard<std::mutex> lock(*_spacialLocks[lockCoords[0]][lockCoords[1]][lockCoords[2]]);
              f->AoSFunctor(p2, p1halo, newton3);
            },
            min, max, IteratorBehavior::owned);
      }
    }
  });
#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerBufferContainer.stop();
  timerPBufferPBuffer.start();
#endif
  // 3. particleBuffer with itself and all other buffers

  // All (halo-)buffer interactions shall happen vectorized, hence, load all buffer data into SoAs
  for (auto &buffer : particleBuffers) {
    f->SoALoader(buffer, buffer._particleSoABuffer, 0);
  }
  for (auto &buffer : haloParticleBuffers) {
    f->SoALoader(buffer, buffer._particleSoABuffer, 0);
  }

#ifdef AUTOPAS_OPENMP
  // Use task parallelism for better dependency resolution than locks
#pragma omp parallel
#pragma omp single
#endif
  {
    for (size_t i = 0; i < particleBuffers.size(); ++i) {
      auto *particleBufferSoA = &particleBuffers[i]._particleSoABuffer;
#ifdef AUTOPAS_OPENMP
#pragma omp task depend(mutexinoutset : particleBufferSoA)
#endif
      f->SoAFunctorSingle(*particleBufferSoA, newton3);
    }

    if constexpr (newton3) {
      for (size_t i = 0; i < particleBuffers.size(); ++i) {
        auto *particleBufferSoAA = &particleBuffers[i]._particleSoABuffer;
        for (size_t j = i + 1; j < particleBuffers.size(); ++j) {
          auto *particleBufferSoAB = &particleBuffers[j]._particleSoABuffer;
#ifdef AUTOPAS_OPENMP
#pragma omp task depend(mutexinoutset : particleBufferSoAA, particleBufferSoAB)
#endif
          f->SoAFunctorPair(*particleBufferSoAA, *particleBufferSoAB, true);
        }
      }
    } else {
      for (size_t i = 0; i < particleBuffers.size(); ++i) {
        auto *particleBufferSoAA = &particleBuffers[i]._particleSoABuffer;
        for (size_t j = 0; j < particleBuffers.size(); ++j) {
          if (i == j) {
            continue;
          }
          auto *particleBufferSoAB = &particleBuffers[j]._particleSoABuffer;
#ifdef AUTOPAS_OPENMP
#pragma omp task depend(mutexinoutset : particleBufferSoAA)
#endif
          f->SoAFunctorPair(*particleBufferSoAA, *particleBufferSoAB, false);
        }
      }
    }
  }

#if SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_TRACE
  timerPBufferPBuffer.stop();
  timerPBufferHBuffer.start();
#endif
// 4. particleBuffer with haloParticleBuffer
#ifdef AUTOPAS_OPENMP
  // Here, phase / color based parallelism turned out to be more efficient than tasks
#pragma omp parallel
#endif
  for (int interactionOffset = 0; interactionOffset < haloParticleBuffers.size(); ++interactionOffset) {
#ifdef AUTOPAS_OPENMP
#pragma omp for
#endif
    for (size_t i = 0; i < particleBuffers.size(); ++i) {
      auto &particleBufferSoA = particleBuffers[i]._particleSoABuffer;
      auto &haloBufferSoA =
          haloParticleBuffers[(i + interactionOffset) % haloParticleBuffers.size()]._particleSoABuffer;
      f->SoAFunctorPair(particleBufferSoA, haloBufferSoA, newton3);
    }
  }
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
template <class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3, bool inTuningPhase>
void AutoTuner<Particle>::iteratePairwiseTemplateHelper(PairwiseFunctor *f, bool doListRebuild,
                                                        std::vector<FullParticleCell<Particle>> &particleBuffer,
                                                        std::vector<FullParticleCell<Particle>> &haloParticleBuffer) {
  auto containerPtr = getContainer();
  AutoPasLog(DEBUG, "Iterating with configuration: {} tuning: {}",
             _tuningStrategy->getCurrentConfiguration().toString(), inTuningPhase ? "true" : "false");

  auto traversal = autopas::utils::withStaticCellType<Particle>(
      containerPtr->getParticleCellTypeEnum(), [&](auto particleCellDummy) {
        return TraversalSelector<decltype(particleCellDummy)>::template generateTraversal<PairwiseFunctor, dataLayout,
                                                                                          useNewton3>(
            _tuningStrategy->getCurrentConfiguration().traversal, *f, containerPtr->getTraversalSelectorInfo());
      });

  if (not traversal->isApplicable()) {
    autopas::utils::ExceptionHandler::exception(
        "Error: Trying to execute a traversal that is not applicable. Two common reasons:\n"
        "1. The search space is trivial, but no traversals are applicable. \n"
        "2. You are using multiple functors and one of the not first is not supporting all configuration options of "
        "the first.\n"
        "Config: {}\n"
        "Current functor: {}",
        _tuningStrategy->getCurrentConfiguration().toString(), typeid(*f).name());
  }

  autopas::utils::Timer timerTotal;
  autopas::utils::Timer timerRebuild;
  autopas::utils::Timer timerIteratePairwise;
  autopas::utils::Timer timerRemainderTraversal;
  timerTotal.start();

  f->initTraversal();
  if (doListRebuild) {
    timerRebuild.start();
    containerPtr->rebuildNeighborLists(traversal.get());
    timerRebuild.stop();
  }
  timerIteratePairwise.start();
  containerPtr->iteratePairwise(traversal.get());
  timerIteratePairwise.stop();

  timerRemainderTraversal.start();
  doRemainderTraversal<useNewton3>(f, containerPtr, particleBuffer, haloParticleBuffer);
  timerRemainderTraversal.stop();
  f->endTraversal(useNewton3);

  timerTotal.stop();

  auto bufferSizeListing = [](const auto &buffers) -> std::string {
    std::stringstream ss;
    size_t sum = 0;
    for (const auto &buffer : buffers) {
      ss << buffer.numParticles() << ", ";
      sum += buffer.numParticles();
    }
    ss << " Total: " << sum;
    return ss.str();
  };
  AutoPasLog(TRACE, "particleBuffer     size : {}", bufferSizeListing(particleBuffer));
  AutoPasLog(TRACE, "haloParticleBuffer size : {}", bufferSizeListing(haloParticleBuffer));
  AutoPasLog(DEBUG, "Container::iteratePairwise took {} ns", timerIteratePairwise.getTotalTime());
  AutoPasLog(DEBUG, "RemainderTraversal         took {} ns", timerRemainderTraversal.getTotalTime());
  AutoPasLog(DEBUG, "RebuildNeighborLists       took {} ns", timerRebuild.getTotalTime());
  AutoPasLog(DEBUG, "AutoTuner::iteratePairwise took {} ns", timerTotal.getTotalTime());

  _iterationLogger.logIteration(getCurrentConfig(), _iteration, inTuningPhase, timerIteratePairwise.getTotalTime(),
                                timerRemainderTraversal.getTotalTime(), timerRebuild.getTotalTime(),
                                timerTotal.getTotalTime());

  // if tuning execute with time measurements
  if (inTuningPhase) {
    if (f->isRelevantForTuning()) {
      addTimeMeasurement(timerTotal.getTotalTime(), doListRebuild);
    } else {
      AutoPasLog(TRACE, "Skipping adding of time measurement because functor is not marked relevant.");
    }
  }
}

template <class Particle>
template <class PairwiseFunctor>
bool AutoTuner<Particle>::tune(PairwiseFunctor &pairwiseFunctor) {
  bool stillTuning = true;

  // need more samples; keep current config
  if (getCurrentNumSamples() < _maxSamples) {
    return stillTuning;
  }
  utils::Timer tuningTimer;
  tuningTimer.start();

  // first tuning iteration -> reset to first config
  if (_iterationsSinceTuning == _tuningInterval) {
    if (auto *mpiStrategy = dynamic_cast<MPIParallelizedStrategy *>(_tuningStrategy.get())) {
      const std::pair<double, double> smoothedHomogeneityAndMaxDensity{
          autopas::OptimumSelector::medianValue(_homogeneitiesOfLastTenIterations),
          autopas::OptimumSelector::medianValue(_maxDensitiesOfLastTenIterations)};
      mpiStrategy->reset<Particle>(_iteration, getContainer(), smoothedHomogeneityAndMaxDensity,
                                   _mpiTuningMaxDifferenceForBucket, _mpiTuningWeightForMaxDensity);
    } else {
      _tuningStrategy->reset(_iteration);
    }
    // only used temporarily of evaluation of smoothing
    if (_tuningStrategy->smoothedHomogeneityAndMaxDensityNeeded()) {
      int rank{0};
#ifdef AUTOPAS_INTERNODE_TUNING
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
      AutoPasLog(DEBUG, "Calculating homogeneities took added up {} ns on rank {}.",
                 _timerCalculateHomogeneity.getTotalTime(), rank);
      _homogeneitiesOfLastTenIterations.erase(_homogeneitiesOfLastTenIterations.begin(),
                                              _homogeneitiesOfLastTenIterations.end());
      _maxDensitiesOfLastTenIterations.erase(_maxDensitiesOfLastTenIterations.begin(),
                                             _maxDensitiesOfLastTenIterations.end());
    }
  } else {  // enough samples -> next config
    stillTuning = _tuningStrategy->tune();
  }

  // repeat as long as traversals are not applicable or we run out of configs
  while (true) {
    auto &currentConfig = getCurrentConfig();
    // check if newton3 works with this functor and remove config if not
    if ((currentConfig.newton3 == Newton3Option::enabled and not pairwiseFunctor.allowsNewton3()) or
        (currentConfig.newton3 == Newton3Option::disabled and not pairwiseFunctor.allowsNonNewton3())) {
      AutoPasLog(WARN, "Configuration with newton 3 {} called with a functor that does not support this!",
                 currentConfig.newton3.to_string());

      _tuningStrategy->removeN3Option(currentConfig.newton3);
    } else {
      if (configApplicable(currentConfig, pairwiseFunctor)) {
        // we found a valid config!
        break;
      } else {
        AutoPasLog(DEBUG, "Skip not applicable configuration {}", currentConfig.toString());
        stillTuning = _tuningStrategy->tune(true);
      }
    }
  }
  // samples should only be cleared if we are still tuning, see `if (_samples.size() < _maxSamples)` from before.
  if (stillTuning) {
    // samples are no longer needed. Delete them here so willRebuild() works as expected.
    _samplesNotRebuildingNeighborLists.clear();
    _samplesRebuildingNeighborLists.clear();
  }
  tuningTimer.stop();
  // when a tuning result is found log it
  if (not stillTuning) {
    AutoPasLog(DEBUG, "Selected Configuration {}", getCurrentConfig().toString());
    _tuningResultLogger.logTuningResult(getCurrentConfig(), _iteration, tuningTimer.getTotalTime());
  }
  AutoPasLog(DEBUG, "Tuning took {} ns.", tuningTimer.getTotalTime());
  _iterationLogger.logTimeTuning(tuningTimer.getTotalTime());

  selectCurrentContainer();
  return stillTuning;
}

template <class Particle>
template <class PairwiseFunctor>
bool AutoTuner<Particle>::configApplicable(const Configuration &conf, PairwiseFunctor &pairwiseFunctor) {
  auto allContainerTraversals = compatibleTraversals::allCompatibleTraversals(conf.container);
  if (allContainerTraversals.find(conf.traversal) == allContainerTraversals.end()) {
    // container and traversal mismatch
    return false;
  }

  _containerSelector.selectContainer(
      conf.container, ContainerSelectorInfo(conf.cellSizeFactor, _verletSkinPerTimestep, _rebuildFrequency,
                                            _verletClusterSize, conf.loadEstimator));
  auto traversalInfo = _containerSelector.getCurrentContainer()->getTraversalSelectorInfo();

  auto containerPtr = getContainer();

  return autopas::utils::withStaticCellType<Particle>(
      containerPtr->getParticleCellTypeEnum(), [&](auto particleCellDummy) {
        return TraversalSelector<decltype(particleCellDummy)>::template generateTraversal<PairwiseFunctor>(
                   conf.traversal, pairwiseFunctor, traversalInfo, conf.dataLayout, conf.newton3)
            ->isApplicable();
      });
}

template <class Particle>
const Configuration &AutoTuner<Particle>::getCurrentConfig() const {
  return _tuningStrategy->getCurrentConfiguration();
}

template <class Particle>
void AutoTuner<Particle>::addTimeMeasurement(long time, bool neighborListRebuilt) {
  const auto &currentConfig = _tuningStrategy->getCurrentConfiguration();
  if (getCurrentNumSamples() < _maxSamples) {
    AutoPasLog(TRACE, "Adding sample.");
    if (neighborListRebuilt) {
      _samplesRebuildingNeighborLists.push_back(time);
    } else {
      _samplesNotRebuildingNeighborLists.push_back(time);
    }
    // if this was the last sample:
    if (getCurrentNumSamples() == _maxSamples) {
      auto &evidenceCurrentConfig = _evidence[currentConfig];

      const long reducedValue = estimateRuntimeFromSamples();

      evidenceCurrentConfig.emplace_back(_iteration, reducedValue);

      // smooth evidence to remove high outliers. If smoothing results in a higher value use the original value.
      const auto smoothedValue = std::min(reducedValue, smoothing::smoothLastPoint(evidenceCurrentConfig, 5));

      // replace collected evidence with smoothed value to improve next smoothing
      evidenceCurrentConfig.back().second = smoothedValue;

      _tuningStrategy->addEvidence(smoothedValue, _iteration);

      // print config, times and reduced value
      if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
        std::ostringstream ss;
        // print config
        ss << currentConfig.toString() << " : ";
        // print all timings
        ss << utils::ArrayUtils::to_string(_samplesRebuildingNeighborLists, " ",
                                           {"With rebuilding neighbor lists [ ", " ]"});
        ss << utils::ArrayUtils::to_string(_samplesNotRebuildingNeighborLists, " ",
                                           {"Without rebuilding neighbor lists [ ", " ]"});
        ss << " Smoothed value: " << smoothedValue;
        AutoPasLog(DEBUG, "Collected times for  {}", ss.str());
      }
      _tuningDataLogger.logTuningData(currentConfig, _samplesRebuildingNeighborLists,
                                      _samplesNotRebuildingNeighborLists, _iteration, reducedValue, smoothedValue);
    }
  }
}
}  // namespace autopas
