/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.2018
 */

#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

#include "autopas/options/TuningMetricOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/selectors/OptimumSelector.h"
#include "autopas/tuning/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/utils/Evidence.h"
#include "autopas/utils/RaplMeter.h"
#include "autopas/utils/Timer.h"
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
 */
class AutoTuner {
 public:
  /**
   * Constructor for the AutoTuner that generates all configurations from the given options.
   * @param tuningStrategy Object implementing the modelling and exploration of a search space.
   * @param MPITuningMaxDifferenceForBucket For MPI-tuning: Maximum of the relative difference in the comparison metric
   * for two ranks which exchange their tuning information.
   * @param MPITuningWeightForMaxDensity For MPI-tuning: Weight for maxDensity in the calculation for bucket
   * distribution.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningMetric Metric used to rate configurations (time or energy).
   * @param tuningInterval Number of time steps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   * @param rebuildFrequency The rebuild frequency this AutoPas instance uses.
   * @param outputSuffix Suffix for all output files produced by this class.
   * @param useTuningStrategyLoggerProxy Whether to use the tuning strategy logger proxy to log tuning information.
   */
  AutoTuner(std::unique_ptr<TuningStrategyInterface> tuningStrategy, double MPITuningMaxDifferenceForBucket,
            double MPITuningWeightForMaxDensity, SelectorStrategyOption selectorStrategy,
            TuningMetricOption tuningMetric, unsigned int tuningInterval, unsigned int maxSamples,
            unsigned int rebuildFrequency, const std::string &outputSuffix = "",
            bool useTuningStrategyLoggerProxy = false);

  /**
   * Move assignment operator
   * @param other
   * @return
   */
  AutoTuner &operator=(AutoTuner &&other) noexcept;

  /**
   * @copydoc autopas::AutoPas::forceRetune()
   */
  void forceRetune();

  /**
   * Getter to a modifiable reference to the current tuning strategy.
   * @return Non-const reference to the underlying tuning strategy.
   */
  [[nodiscard]] TuningStrategyInterface &getTuningStrategy() const { return *_tuningStrategy; }

  /**
   * Getter for the primary metric for tuning.
   * @return
   */
  const TuningMetricOption &getTuningMetric() const { return _tuningMetric; }

  /**
   * Determines what live infos are needed and resets the strategy upon the start of a new tuning phase.
   *
   * @note The live info is not gathered here because then we would need
   *
   * @return Tuple indicating what is needed before the next call to tune: tuple<liveInfo, homogeneity>
   */
  std::tuple<bool, bool> prepareIteration() {
    // Flag if this is the first iteration in a new tuning phase
    const bool startOfTuningPhase = _iterationsSinceTuning == _tuningInterval;

    // first tuning iteration -> reset everything
    if (startOfTuningPhase) {
      // call the appropriate version of reset
      if (auto *mpiStrategy = dynamic_cast<MPIParallelizedStrategy *>(_tuningStrategy.get())) {
        const std::pair<double, double> smoothedHomogeneityAndMaxDensity{
            autopas::OptimumSelector::medianValue(_homogeneitiesOfLastTenIterations),
            autopas::OptimumSelector::medianValue(_maxDensitiesOfLastTenIterations)};
        mpiStrategy->reset(_iteration, smoothedHomogeneityAndMaxDensity, _mpiTuningMaxDifferenceForBucket,
                           _mpiTuningWeightForMaxDensity);
      } else {
        _tuningStrategy->reset(_iteration);
      }

      // Homogeneity was calculated directly before the tuning phase so reset it now.
      if (_tuningStrategy->smoothedHomogeneityAndMaxDensityNeeded()) {
        AutoPasLog(DEBUG, "Calculating homogeneities over 10 iterations took in total {} ns on rank {}.",
                   _timerCalculateHomogeneity.getTotalTime(), []() {
                     int rank{0};
                     AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
                     return rank;
                   });
        _homogeneitiesOfLastTenIterations.clear();
        _maxDensitiesOfLastTenIterations.clear();
      }
    }

    // if necessary, we need to collect live info in the first tuning iteration
    const bool needsLiveInfo = startOfTuningPhase and _tuningStrategy->needsLiveInfo();

    // calc homogeneity if needed, and we are within 10 iterations of the next tuning phase
    const size_t numIterationsForHomogeneity = 10;
    const bool needsHomogeneityAndMaxDensity =
        _tuningStrategy->smoothedHomogeneityAndMaxDensityNeeded() and
        _iterationsSinceTuning > _tuningInterval - numIterationsForHomogeneity and
        _iterationsSinceTuning <= _tuningInterval;

    return {needsLiveInfo, needsHomogeneityAndMaxDensity};
  }

  /**
   * Increase internal iteration counters by one.
   */
  void bumpIterationCounters();

  /**
   * Returns whether rebuildNeighborLists() will be triggered in the next call to iteratePairwise().
   * This might also indicate a container change.
   *
   * @return True if the the current iteration counters indicate a rebuild in the next iteration.
   */
  bool willRebuildNeighborLists() const;

  /**
   * Get the currently selected configuration.
   * @return
   */
  [[nodiscard]] const Configuration &getCurrentConfig() const;

  /**
   * Ask the tuner for the next configuration to use.
   * This either returns the already selected config or triggers a step of the tuning process.
   * @return Tuple<Next configuration to use, still tuning>.
   */
  [[nodiscard]] std::tuple<Configuration, bool> getNextConfig();

  /**
   * Tell the tuner that the given config is not applicable.
   *
   * @note The applicability checking logic was moved out of the tuner because it needed the container,
   * thus raising the compile complexity.
   *
   * @tparam PairwiseFunctor
   * @param rejectedConfig
   * @param indefinitely Whether the given config should be completely removed from the search space (aka rejected
   * indefinitely).
   * @return
   */
  [[nodiscard]] std::tuple<Configuration, bool> rejectConfig(const Configuration &rejectedConfig, bool indefinitely) {
    if (searchSpaceIsTrivial()) {
      utils::ExceptionHandler::exception("Rejected the only configuration in the search space!\n{}",
                                         rejectedConfig.toString());
    }

    if (indefinitely) {
      // delete rejected config from the search space and assemble new sets.
      _searchSpace.erase(rejectedConfig);
    }

    const auto stillTuning = _tuningStrategy->tune(true);
    return {getCurrentConfig(), stillTuning};
  }

  /**
   * Initialize the container specified by the TuningStrategy.
   */
  bool searchSpaceIsTrivial();

  /**
   * Log the collected data and if we are at the end of a tuning phase the result to files.
   * @param conf
   * @param tuningIteration
   * @param tuningTime
   */
  void logIteration(const Configuration &conf, bool tuningIteration, long tuningTime);

  /**
   * Initialize rapl meter.
   * @return True if energy measurements are possible on this system.
   */
  bool initEnergy();

  /**
   * Reset the rapl meter to prepare for a new measurement.
   * @return True if energy measurements are possible on this system.
   */
  bool resetEnergy();

  /**
   * Take an energy measurement.
   * @return Tuple<PsysEnergy, PkgEnergy, RamEnergy, TotalEnergy>
   */
  std::tuple<double, double, double, long> sampleEnergy();

  /**
   * Save the runtime of a given traversal.
   *
   * Samples are collected and reduced to one single value according to _selectorStrategy. Only then the value is passed
   * on to the tuning strategy. This function expects that samples of the same configuration are taken consecutively.
   * The sample argument is a long because std::chrono::duration::count returns a long.
   *
   * @param sample
   * @param neighborListRebuilt If the neighbor list as been rebuilt during the given time.
   */
  void addMeasurement(long sample, bool neighborListRebuilt);

  /**
   * Adds measurements of homogeneity and maximal density to the vector of measurements.
   * @param homogeneity
   * @param maxDensity
   * @param time Time it took to obtain these measurements.
   */
  void addHomogeneityAndMaxDensity(double homogeneity, double maxDensity, long time);

 private:
  /**
   * Measures consumed energy for tuning
   */
  utils::RaplMeter _raplMeter;

  /**
   * Total number of collected samples. This is the sum of the sizes of all sample vectors.
   * @return Sum of sizes of sample vectors.
   */
  size_t getCurrentNumSamples() const;

  /**
   * Estimate the runtime from the current samples according to the SelectorStrategy and rebuild frequency.
   * Samples are weighted so that we normalize to the expected number of (non-)rebuild iterations and then divide by the
   * rebuild frequency.
   * @return estimate time for one iteration
   */
  [[nodiscard]] long estimateRuntimeFromSamples() const;

  /**
   * Tune available algorithm configurations.
   *
   * It is assumed this function is only called for relevant functors and that at least two configurations are allowed.
   * When in tuning phase selects next config to test. At the end of the tuning phase select optimum.
   * The function returns true if the selected config is not yet the optimum but something that should be sampled.
   *
   * @return true iff still in tuning phase.
   */
  bool tune();

  SelectorStrategyOption _selectorStrategy;
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  /**
   * Counter for the current simulation iteration.
   */
  size_t _iteration;

  /**
   * Number of iterations between two tuning phases.
   */
  size_t _tuningInterval;

  /**
   * Number of iterations since the end of the last tuning phase.
   */
  size_t _iterationsSinceTuning;

  /**
   * Metric to use for tuning.
   */
  TuningMetricOption _tuningMetric;

  /**
   * Is energy measurement possible.
   * Initialize as true and check in the constructor if it is indeed possible.
   */
  bool _energyMeasurementPossible;

  /**
   * The rebuild frequency this instance of AutoPas uses.
   */
  unsigned int _rebuildFrequency;

  /**
   * How many times each configuration should be tested.
   */
  size_t _maxSamples;

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
  std::vector<double> _homogeneitiesOfLastTenIterations{};

  /**
   * Buffer for the homogeneities of the last ten Iterations
   */
  std::vector<double> _maxDensitiesOfLastTenIterations{};

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
  std::map<Configuration, std::vector<Evidence>> _searchSpace;

  /**
   * Timer used to determine how much time is wasted by calculating multiple homogeneities for smoothing
   * Only used temporarily
   */
  autopas::utils::Timer _timerCalculateHomogeneity;

  TuningResultLogger _tuningResultLogger;
  TuningDataLogger _tuningDataLogger;
};
}  // namespace autopas