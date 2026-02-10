/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.2018
 */

#pragma once

#include <cstddef>
#include <memory>
#include <set>
#include <tuple>

#include "autopas/options/TuningMetricOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/LiveInfo.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/utils/AutoTunerInfo.h"
#include "autopas/utils/EnergySensor.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/TuningDataLogger.h"
#include "autopas/utils/logging/TuningResultLogger.h"

namespace autopas {

/**
 * This class manages all logic related to the auto tuning mechanic. This involves:
 *  - Managing the search space and configQueue.
 *  - Managing and triggering all active tuning strategies.
 *  - Managing all collected evidence and passing it to the strategies.
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
   * Type for the member holding all tuning strategies.
   */
  using TuningStrategiesListType = std::vector<std::unique_ptr<TuningStrategyInterface>>;
  /**
   * Type for the search space holding all possible configurations.
   */
  using SearchSpaceType = std::set<Configuration>;

  /**
   * Constructor for the AutoTuner that generates all configurations from the given options.
   * @param tuningStrategies Vector of object implementing the modelling and exploration of a search space. Will be
   * moved into the tuner.
   * @param searchSpace All possible configurations.
   * @param autoTunerInfo Struct containing more configuration information.
   * @param rebuildFrequency The number of iterations after which the neighbor lists are rebuilt.
   * @param outputSuffix Suffix for all output files produced by this object.
   */
  AutoTuner(TuningStrategiesListType &tuningStrategies, const SearchSpaceType &searchSpace,
            const AutoTunerInfo &autoTunerInfo, unsigned int rebuildFrequency, const std::string &outputSuffix);

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
   * Getter for the primary metric for tuning.
   * @return
   */
  const TuningMetricOption &getTuningMetric() const;

  /**
   * Pass live info on to all tuning strategies.
   * @param liveInfo
   */
  void receiveLiveInfo(const LiveInfo &liveInfo);

  /**
   * Returns true if the AutoTuner is about to calculate the first interactions of a tuning phase (i.e. the first
   * iteration), before tuneConfiguration() has been called.
   * @return
   */
  bool isStartOfTuningPhase() const;

  /**
   * Returns true if the AutoTuner is within 10 iterations of the start of a tuning phase.
   * @return
   */
  bool tuningPhaseAboutToBegin() const;

  /**
   * Returns true if the AutoTuner needs live info. This occurs if any strategy requires this and AutoPas is beginning
   * a tuning phase or if a strategy requires domain similarity statistics (taken from LiveInfo) and AutoPas is within
   * 10 iterations of a tuning phase.
   * @return True if the AutoTuner needs live info.
   */
  [[nodiscard]] bool needsLiveInfo() const;

  /**
   * Increase internal iteration counters by one. Should be called at the end of an iteration.
   * @param needToWait If tuner should wait for other tuners.
   */
  void bumpIterationCounters(bool needToWait = false);

  /**
   * Returns whether rebuildNeighborLists() should be triggered in the next iteration.
   * This indicates a configuration change.
   * In the non-tuning phase, the rebuildNeighborLists() is triggered in LogicHandler.
   * @return True if the current iteration counters indicate a rebuild in the next iteration due to a configuration
   * change.
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
   * Since this operation might change the suggestion what configuration to try next, this next suggestion is returned.
   *
   * @note The applicability checking logic was moved out of the tuner because it needed the container,
   * thus raising the compile complexity.
   *
   * @param rejectedConfig
   * @param indefinitely Whether the given config should be completely removed from the search space (aka rejected
   * indefinitely).
   * @return Tuple<Next configuration to use, still tuning>.
   */
  [[nodiscard]] std::tuple<Configuration, bool> rejectConfig(const Configuration &rejectedConfig, bool indefinitely);

  /**
   * Indicator function whether the search space consists of exactly one configuration.
   * @return
   */
  bool searchSpaceIsTrivial() const;

  /**
   * Indicator function whether the search space has no configurations in it.
   * @return
   */
  bool searchSpaceIsEmpty() const;

  /**
   * After a tuning phase has finished, write the result to a file.
   * @param tuningIteration
   * @param tuningTime
   */
  void logTuningResult(bool tuningIteration, long tuningTime) const;

  /**
   * Initialize pmt sensor.
   * @return True if energy measurements are enabled and possible.
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
   * Adds domain similarity statistics to a vector of measurements, which can be smoothed for use in MPI Tuning to find
   * similar domains.
   * @param pdBinDensityStdDev particle-dependent bin density standard deviation. See LiveInfo::gather for more
   * information.
   * @param pdBinMaxDensity particle-dependent bin maximum density. See LiveInfo::gather for more information.
   */
  void addDomainSimilarityStatistics(double pdBinDensityStdDev, double pdBinMaxDensity);

  /**
   * Getter for the current queue of configurations.
   * @return
   */
  const std::vector<Configuration> &getConfigQueue() const;

  /**
   * Get the list of tuning strategies that are used.
   * @return
   */
  const std::vector<std::unique_ptr<TuningStrategyInterface>> &getTuningStrategies() const;

  /**
   * Indicate if the tuner considers itself currently in a tuning phase according to its internal counters.
   * @return
   */
  bool inTuningPhase() const;

  /**
   * Indicate if the tuner is in the first iteration of a tuning phase.
   * @return
   */
  bool inFirstTuningIteration() const;

  /**
   * Indicates if the tuner is in the last iteration of the tuning phase.
   * @return
   */
  bool inLastTuningIteration() const;

  /**
   * Indicates whether the tuner is in the iteration corresponding to the last sample of the first configuration in the
   * current tuning phase.
   * @return
   */
  bool inFirstConfigurationLastSample() const;

  /**
   * Getter for the internal evidence collection.
   * @return
   */
  const EvidenceCollection &getEvidenceCollection() const;

  /**
   * Returns whether the AutoTuner can take energy measurements.
   * @return
   */
  bool canMeasureEnergy() const;

  /**
   * Sets the _rebuildFrequency. This is the average number of iterations per rebuild.
   * This is used to dynamically change the _rebuildFrequency based on estimate in case of dynamic containers.
   * @param rebuildFrequency Current rebuild frequency in this instance of autopas, used by autopas for weighing rebuild
   * and non-rebuild iterations
   */
  void setRebuildFrequency(double rebuildFrequency);

  /**
   * Checks whether the current configuration performs so poorly that it shouldn't be resampled further within this
   * tuning phase. If the currently sampled configuration is worse than the current best configuration by more than the
   * earlyStoppingFactor factor, it will not be sampled again this tuning phase. Uses the _estimateRuntimeFromSamples()
   * function to estimate the runtimes.
   */
  void checkEarlyStoppingCondition();

 private:
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
  bool tuneConfiguration();

  /**
   * If it is the end of the tuning phase, determine the optimal configuration and set this as the configuration to be
   * used until the next tuning phase, as well as setting other relevant class members (_endOfTuningPhase, _isTuning,
   * _samplesRebuildingNeighborLists, _iterationBaseline)
   */
  void handleEndOfTuningPhaseIfRelevant();

  /**
   * Strategy how to reduce the sampled values to one value.
   */
  SelectorStrategyOption _selectorStrategy;

  /**
   * Vector holding all tuning strategies that this auto tuner applies.
   * The strategies are always applied in the order they are in this vector.
   */
  std::vector<std::unique_ptr<TuningStrategyInterface>> _tuningStrategies;

  /**
   * Counter for the current simulation iteration.
   * The first iteration has number 0.
   */
  size_t _iteration{0};

  /**
   * The number of the current tuning phase.
   * If we are currently between phases this is the number of the last phase.
   * The first tuning phase has number 0.
   * See bumpIterationCounters() for more details.
   */
  size_t _tuningPhase{0};

  /**
   * Fixed interval at which tuning phases are started.
   * A tuning phase always starts when _iteration % _tuningInterval == 0.
   */
  size_t _tuningInterval;

  /**
   * Metric to use for tuning.
   */
  TuningMetricOption _tuningMetric;

  /**
   * Flag for whether LOESS Smoothening is used to smoothen the tuning results.
   */
  bool _useLOESSSmoothening;

  /**
   * Is energy measurement possible.
   * Initialize as true and check in the constructor if it is indeed possible.
   */
  bool _energyMeasurementPossible;

  /**
   * The rebuild frequency this instance of AutoPas uses.
   * In the case of dynamic containers, this reflects the current estimate of the average rebuild frequency.
   * As the estimate is an average of integers, it can be a fraction and so double is used here.
   * In static containers, this is set to user-defined rebuild frequency.
   */
  double _rebuildFrequency;

  /**
   * How many times each configuration should be tested.
   */
  size_t _maxSamples;

  /**
   * EarlyStoppingFactor for the auto-tuner. A configuration performing worse than the currently best configuration
   * by more than this factor will not be sampled again this tuning phase.
   */
  double _earlyStoppingFactor;

  /**
   * Flag indicating if any tuning strategy needs the smoothed homogeneity and max density collected.
   */
  bool _needsDomainSimilarityStatistics;

  /**
   * Flag indicating if any tuning strategy needs live infos collected.
   */
  bool _needsLiveInfo;

  /**
   * Buffer for the homogeneities of the last ten Iterations
   */
  std::vector<double> _pdBinDensityStdDevOfLastTenIterations{};

  /**
   * Buffer for the maximum densities of the last ten Iterations.
   */
  std::vector<double> _pdBinMaxDensityOfLastTenIterations{};

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
   * Database of all evidence collected so far.
   */
  EvidenceCollection _evidenceCollection{};

  /**
   * The search space for this tuner.
   */
  SearchSpaceType _searchSpace;

  /**
   * Sorted queue of configurations that should be looked at in this tuning phase.
   * Initially this is the full search space. Tuning strategies then can filter and resort this.
   *
   * @note The next configuration to use is at the END of the vector.
   */
  std::vector<Configuration> _configQueue;

  /**
   * CSV logger for configurations selected at the end of each tuning phase.
   */
  TuningResultLogger _tuningResultLogger;

  /**
   * CSV logger for all samples collected during a tuning phase.
   */
  TuningDataLogger _tuningDataLogger;

  /**
   * Sensor for energy measurement
   */
  utils::EnergySensor _energySensor;

  /**
   * Is set to true during a tuning phase.
   */
  bool _isTuning{false};

  /**
   * Is only set to true for the last iteration of a tuning phase.
   * Specifically, from when tuneConfiguration() selects the optimum until it is reset to false in
   * bumpIterationCounters().
   */
  bool _endOfTuningPhase{false};

  /**
   * Is set to true in forceRetune() to signal a new tuning phase should start outside the regular tuningInterval. Is
   * set back to false in tuneConfiguration()
   */
  bool _forceRetune{false};

  /**
   * Is set to true if the current configuration should not be used to collect further samples because it is
   * significantly slower than the fastest configuration by more than _earlyStoppingFactor.
   */
  bool _earlyStoppingOfResampling{false};

  /**
   * Used only for triggering rebuilds when configurations switch during tuning phases, which occurs when
   * _iterationBaseline % _maxSamples == 0. _iterationBaseline may therefore be modified to "skip" iterations e.g. when
   * early stopping is used."
   */
  size_t _iterationBaseline{0};
};
}  // namespace autopas
