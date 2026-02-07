/**
 * @file AutoTuner.cpp
 * @author F. Gratl
 * @date 21.06.23
 */

#include "AutoTuner.h"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>

#include "autopas/tuning/selectors/OptimumSelector.h"
#include "autopas/tuning/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/tuning/utils/Smoothing.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {
AutoTuner::AutoTuner(TuningStrategiesListType &tuningStrategies, const SearchSpaceType &searchSpace,
                     const AutoTunerInfo &autoTunerInfo, unsigned int rebuildFrequency, const std::string &outputSuffix)
    : _selectorStrategy(autoTunerInfo.selectorStrategy),
      _tuningStrategies(std::move(tuningStrategies)),
      _tuningInterval(autoTunerInfo.tuningInterval),
      _tuningMetric(autoTunerInfo.tuningMetric),
      _useLOESSSmoothening(autoTunerInfo.useLOESSSmoothening),
      _energyMeasurementPossible(initEnergy()),
      _rebuildFrequency(rebuildFrequency),
      _maxSamples(autoTunerInfo.maxSamples),
      _earlyStoppingFactor(autoTunerInfo.earlyStoppingFactor),
      _needsDomainSimilarityStatistics(
          std::transform_reduce(_tuningStrategies.begin(), _tuningStrategies.end(), false, std::logical_or(),
                                [](auto &tuningStrat) { return tuningStrat->needsDomainSimilarityStatistics(); })),
      _needsLiveInfo(std::transform_reduce(_tuningStrategies.begin(), _tuningStrategies.end(), false, std::logical_or(),
                                           [](auto &tuningStrat) { return tuningStrat->needsLiveInfo(); })),
      _samplesNotRebuildingNeighborLists(autoTunerInfo.maxSamples),
      _searchSpace(searchSpace),
      _configQueue(searchSpace.begin(), searchSpace.end()),
      _tuningResultLogger(outputSuffix),
      _tuningDataLogger(autoTunerInfo.maxSamples, outputSuffix),
      _energySensor(utils::EnergySensor(autoTunerInfo.energySensor)) {
  _samplesRebuildingNeighborLists.reserve(autoTunerInfo.maxSamples);
  _pdBinDensityStdDevOfLastTenIterations.reserve(10);
  _pdBinMaxDensityOfLastTenIterations.reserve(10);
  if (_searchSpace.empty()) {
    utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
  }
  // AutoTuner will start out with tuning if there are more than one configuration.
  if (_searchSpace.size() > 1) {
    _isTuning = true;
  }
  AutoPasLog(DEBUG, "Points in search space: {}", _searchSpace.size());
  AutoPasLog(DEBUG, "AutoTuner constructed with LOESS Smoothening {}.", _useLOESSSmoothening ? "enabled" : "disabled");
}

AutoTuner &AutoTuner::operator=(AutoTuner &&other) noexcept {
  _tuningStrategies = std::move(other._tuningStrategies);
  return *this;
}

void AutoTuner::addDomainSimilarityStatistics(double pdBinDensityStdDev, double pdBinMaxDensity) {
  _pdBinDensityStdDevOfLastTenIterations.push_back(pdBinDensityStdDev);
  _pdBinMaxDensityOfLastTenIterations.push_back(pdBinMaxDensity);
}

void AutoTuner::logTuningResult(long tuningTime) const {
  // only log if we are at the end of a tuning phase
  if (_endOfTuningPhase) {
    // This string is part of several older scripts, hence it is not recommended to change it.
    const auto [conf, optimalEvidence] = _evidenceCollection.getLatestOptimalConfiguration();
    AutoPasLog(DEBUG, "Selected Configuration {}", conf.toString());
    _tuningResultLogger.logTuningResult(conf, _iteration, tuningTime, optimalEvidence.value);
  }
}

bool AutoTuner::searchSpaceIsTrivial() const { return _searchSpace.size() == 1; }

bool AutoTuner::searchSpaceIsEmpty() const { return _searchSpace.empty(); }

void AutoTuner::forceRetune() {
  if (inTuningPhase()) {
    AutoPasLog(WARN, "Warning: Currently running tuning phase is aborted a new one is started!");
  }
  _samplesNotRebuildingNeighborLists.resize(_maxSamples);
  _forceRetune = true;
  _requiresRebuild = true;
}

bool AutoTuner::tuneConfiguration() {
  utils::Timer tuningTimer;
  tuningTimer.start();

  _requiresRebuild = false;
  // If we are not (yet) tuning or there is nothing to tune return immediately.
  if (not inTuningPhase()) {
    return false;
  }
  if (getCurrentNumSamples() < _maxSamples and not _earlyStoppingOfResampling) {
    // If we are still collecting samples from one config return immediately.
    return true;
  }

  // We are looking for a new config, so we will require a rebuild
  _requiresRebuild = true;

  // We finished collection samples for this config so remove it from the queue
  _configQueue.pop_back();

  // We plan to test a new config so clear all samples.
  _samplesNotRebuildingNeighborLists.clear();
  _samplesRebuildingNeighborLists.clear();

  // Helper function to reset the ConfigQueue if something wipes it.
  auto restoreConfigQueueIfEmpty = [&](const auto &configQueueBackup, const TuningStrategyOption &stratOpt) {
    if (_configQueue.empty()) {
      _configQueue = configQueueBackup;
      AutoPasLog(WARN, "ConfigQueue wipe by {} detected! Resetting to previous state: (Size={}) {}",
                 stratOpt.to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
    }
  };

  // Determine where in a tuning phase we are
  // If _iterationsInMostRecentTuningPhase >= _tuningInterval the current tuning phase takes more iterations than the
  // tuning interval -> continue tuning
  if (isStartOfTuningPhase()) {
    // CASE: Start of a new tuning phase
    _isTuning = true;
    _forceRetune = false;
    // in the first iteration of a tuning phase we reset all strategies
    // and refill the queue with the complete search space.
    // Reverse the order, because _configQueue is FiLo, and we aim to keep the order for legacy reasons.
    _configQueue.clear();
    _configQueue.reserve(_searchSpace.size());
    std::copy(_searchSpace.rbegin(), _searchSpace.rend(), std::back_inserter(_configQueue));
    AutoPasLog(DEBUG, "ConfigQueue at tuneConfiguration before reset: (Size={}) {}", _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
    // then let the strategies filter and sort it
    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
      const auto configQueueBackup = _configQueue;
      const auto intentionalWipe = tuningStrategy->reset(_iteration, _tuningPhase, _configQueue, _evidenceCollection);
      AutoPasLog(DEBUG, "ConfigQueue after applying {}::reset(): (Size={}) {}",
                 tuningStrategy->getOptionType().to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
      if (not intentionalWipe) {
        restoreConfigQueueIfEmpty(configQueueBackup, tuningStrategy->getOptionType());
      }
    });
  } else {
    // CASE: somewhere in a tuning phase
    _isTuning = true;

    AutoPasLog(DEBUG, "ConfigQueue at tuneConfiguration before optimizeSuggestions: (Size={}) {}", _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
      const auto configQueueBackup = _configQueue;
      const auto intentionalWipe = tuningStrategy->optimizeSuggestions(_configQueue, _evidenceCollection);
      AutoPasLog(DEBUG, "ConfigQueue after applying {}::optimizeSuggestions(): (Size={}) {}",
                 tuningStrategy->getOptionType().to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
      if (not intentionalWipe) {
        restoreConfigQueueIfEmpty(configQueueBackup, tuningStrategy->getOptionType());
      }
    });
  }
  // Apply the Container Constraint (Filter).
  if (_containerConstraint.has_value()) {
    // Remove any configuration that does not match the forced container.
    std::erase_if(_configQueue,
                  [&](const Configuration &conf) { return conf.container != _containerConstraint.value(); });
  }

  // CASE: End of a tuning phase. This is not exclusive to the other cases!
  if (_configQueue.empty()) {
    // If the queue is empty we are done tuning.
    _endOfTuningPhase = true;
    selectBestConfiguration();
  }
  tuningTimer.stop();

  _earlyStoppingOfResampling = false;
  return _isTuning;
}

void AutoTuner::selectBestConfiguration() {
  // If there are still queued up configurations, clear them.
  _configQueue.clear();
  // Find and push_back the optimal configuration for the current container.
  const auto [optConf, optEvidence] = _evidenceCollection.getOptimalConfiguration(_tuningPhase, _containerConstraint);
  _configQueue.push_back(optConf);
  _isTuning = false;
  // Fill up sample buffer to indicate we are not collecting samples anymore.
  _samplesRebuildingNeighborLists.resize(_maxSamples);
}

const Configuration &AutoTuner::getCurrentConfig() const {
  if (_configQueue.empty()) {
    utils::ExceptionHandler::exception(
        "AutoTuner::getCurrentConfig(): Cannot get the current Configuration as the config queue is empty");
  }
  return _configQueue.back();
}

std::tuple<Configuration, bool> AutoTuner::rejectConfig(const Configuration &rejectedConfig, bool indefinitely) {
  if (searchSpaceIsTrivial()) {
    utils::ExceptionHandler::exception("Rejected the only configuration in the search space!\n{}",
                                       rejectedConfig.toString());
  }

  // remove the config from the queue
  std::erase_if(_configQueue, [&](const auto &conf) { return conf == rejectedConfig; });

  if (indefinitely) {
    // delete rejected config from the search space and notify tuning strategies.
    _searchSpace.erase(rejectedConfig);
  }
  std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
    tuningStrategy->rejectConfiguration(rejectedConfig, indefinitely);
    AutoPasLog(DEBUG, "ConfigQueue after applying {}::rejectConfiguration(): (Size={}) {}",
               tuningStrategy->getOptionType().to_string(), _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
  });

  const auto stillTuning = tuneConfiguration();
  return {getCurrentConfig(), stillTuning};
}

void AutoTuner::addMeasurement(long sample, bool neighborListRebuilt) {
  const auto &currentConfig = _configQueue.back();
  // sanity check
  if (getCurrentNumSamples() >= _maxSamples) {
    utils::ExceptionHandler::exception(
        "AutoTuner::addMeasurement(): Trying to add a new measurement to the AutoTuner but there are already enough "
        "for this configuration!\n"
        "tuneConfiguration() should have been called before to process and flush samples.");
  }
  AutoPasLog(TRACE, "Adding sample {} to configuration {}.", sample, currentConfig.toShortString());
  if (neighborListRebuilt) {
    _samplesRebuildingNeighborLists.push_back(sample);
  } else {
    _samplesNotRebuildingNeighborLists.push_back(sample);
  }

  checkEarlyStoppingCondition();

  // if this was the last sample for this configuration:
  //  - calculate the evidence from the collected samples
  //  - log what was collected
  //  - remove the configuration from the queue
  if (getCurrentNumSamples() == _maxSamples or _earlyStoppingOfResampling) {
    const long reducedValue = estimateRuntimeFromSamples();
    _evidenceCollection.addEvidence(currentConfig, {_iteration, _tuningPhase, reducedValue});

    // If LOESS-based smoothening is enabled, use it to smooth evidence to remove high outliers. If smoothing results in
    // a higher value or if LOESS-based smoothening is disabled, use the original value.
    const auto smoothedValue =
        _useLOESSSmoothening
            ? std::min(reducedValue, smoothing::smoothLastPoint(*_evidenceCollection.getEvidence(currentConfig), 5))
            : reducedValue;

    // replace collected evidence with smoothed value to improve next smoothing
    _evidenceCollection.modifyLastEvidence(currentConfig).value = smoothedValue;

    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
      tuningStrategy->addEvidence(getCurrentConfig(), _evidenceCollection.modifyLastEvidence(currentConfig));
    });

    // print config, times and reduced value
    AutoPasLog(
        DEBUG, "Collected {} for {}",
        [&]() {
          switch (this->_tuningMetric) {
            case TuningMetricOption::time:
              return "times";
            case TuningMetricOption::energy:
              return "energy consumption";
          }
          autopas::utils::ExceptionHandler::exception("AutoTuner::addMeasurement(): Unknown tuning metric.");
          return "Unknown tuning metric";
        }(),
        [&]() {
          std::ostringstream ss;
          // print config
          ss << currentConfig << " : ";
          // print all timings
          ss << utils::ArrayUtils::to_string(_samplesRebuildingNeighborLists, " ",
                                             {"With rebuilding neighbor lists [ ", " ] "});
          ss << utils::ArrayUtils::to_string(_samplesNotRebuildingNeighborLists, " ",
                                             {"Without rebuilding neighbor lists [ ", " ] "});
          ss << "Smoothed value: " << smoothedValue;
          return ss.str();
        }());

    auto samplesRebuildingNeighborLists = _samplesRebuildingNeighborLists;
    auto samplesNotRebuildingNeighborLists = _samplesNotRebuildingNeighborLists;

    if (_earlyStoppingOfResampling) {
      // pad sample vectors to length of maxSamples to ensure correct logging
      samplesNotRebuildingNeighborLists.resize(_maxSamples - samplesRebuildingNeighborLists.size(), -1);
    }

    _tuningDataLogger.logTuningData(currentConfig, samplesRebuildingNeighborLists, samplesNotRebuildingNeighborLists,
                                    _iteration, reducedValue, smoothedValue, _rebuildFrequency);
  }
}

void AutoTuner::bumpIterationCounters() {
  ++_iteration;
  AutoPasLog(DEBUG, "Iteration: {}", _iteration);
  _endOfTuningPhase = false;
}

bool AutoTuner::willRebuildNeighborLists() const { return _requiresRebuild; }

bool AutoTuner::initEnergy() {
  // Try to initialize the raplMeter
  return _energySensor.init(_tuningMetric == TuningMetricOption::energy);
}

bool AutoTuner::resetEnergy() { return _energySensor.startMeasurement(); }

std::tuple<double, double, double, long> AutoTuner::sampleEnergy() {
  _energySensor.endMeasurement();
  return {_energySensor.getWatts(), _energySensor.getJoules(), _energySensor.getEnergyDeltaT(),
          _energySensor.getNanoJoules()};
}

size_t AutoTuner::getCurrentNumSamples() const {
  return _samplesNotRebuildingNeighborLists.size() + _samplesRebuildingNeighborLists.size();
}

long AutoTuner::estimateRuntimeFromSamples() const {
  // reduce samples for rebuild and non-rebuild iterations with the given selector strategy
  const auto reducedValueBuilding =
      autopas::OptimumSelector::optimumValue(_samplesRebuildingNeighborLists, _selectorStrategy);
  // if there is no data for the non rebuild iterations we have to assume them taking the same time as rebuilding ones
  // this might neither be a good estimate nor fair but the best we can do
  const auto reducedValueNotBuilding =
      _samplesNotRebuildingNeighborLists.empty()
          ? reducedValueBuilding
          : autopas::OptimumSelector::optimumValue(_samplesNotRebuildingNeighborLists, _selectorStrategy);

  // Calculate weighted average as if there was exactly one sample for each iteration in the rebuild interval.
  return (reducedValueBuilding + (_rebuildFrequency - 1) * reducedValueNotBuilding) / _rebuildFrequency;
}

bool AutoTuner::isStartOfTuningPhase() const {
  return (_iteration % _tuningInterval == 0 and not _isTuning) or _forceRetune;
}

bool AutoTuner::tuningPhaseAboutToBegin() const { return _iteration % _tuningInterval > _tuningInterval - 10; }

bool AutoTuner::needsLiveInfo() const {
  return (isStartOfTuningPhase() and _needsLiveInfo) or
         (tuningPhaseAboutToBegin() and _needsDomainSimilarityStatistics);
}

const std::vector<Configuration> &AutoTuner::getConfigQueue() const { return _configQueue; }

const std::vector<std::unique_ptr<TuningStrategyInterface>> &AutoTuner::getTuningStrategies() const {
  return _tuningStrategies;
}

void AutoTuner::receiveLiveInfo(const LiveInfo &liveInfo) {
  // Handle Live Info processing before the tuning phase
  if (_needsDomainSimilarityStatistics and tuningPhaseAboutToBegin()) {
    const auto particleDependentBinDensityStdDev = liveInfo.get<double>("particleDependentBinDensityStdDev");
    const auto particleDependentBinMaxDensity = liveInfo.get<double>("particleDependentBinMaxDensity");
    addDomainSimilarityStatistics(particleDependentBinDensityStdDev, particleDependentBinMaxDensity);
  }
  // Handling at start of tuning phase
  if ((_needsLiveInfo or _needsDomainSimilarityStatistics) and isStartOfTuningPhase()) {
    for (auto &tuningStrategy : _tuningStrategies) {
      if (tuningStrategy->needsLiveInfo()) {
        tuningStrategy->receiveLiveInfo(liveInfo);
      }
      if (tuningStrategy->needsDomainSimilarityStatistics()) {
        const auto pdBinDensityStdDevSmoothed = OptimumSelector::medianValue(_pdBinDensityStdDevOfLastTenIterations);
        const auto pdBinMaxDensitySmoothed = OptimumSelector::medianValue(_pdBinMaxDensityOfLastTenIterations);
        tuningStrategy->receiveDomainSimilarityStatistics(pdBinDensityStdDevSmoothed, pdBinMaxDensitySmoothed);
      }
    }
  }
  if (_needsDomainSimilarityStatistics and _forceRetune and isStartOfTuningPhase()) {
    // If we have sent statistics because of a forced retune, throw a warning
    AutoPasLog(WARN,
               "Due to a forced retune, domain similarity statistics, used by at least one tuning strategy, may"
               "not be correct.");
  }
}

const TuningMetricOption &AutoTuner::getTuningMetric() const { return _tuningMetric; }

bool AutoTuner::inTuningPhase() const { return (_isTuning or _forceRetune) and not searchSpaceIsTrivial(); }

bool AutoTuner::inFirstTuningIteration() const { return (_iteration % _tuningInterval == 0); }

bool AutoTuner::inLastTuningIteration() const { return _endOfTuningPhase; }

bool AutoTuner::inFirstConfigurationLastSample() const { return (_iteration % _tuningInterval == _maxSamples - 1); }

const EvidenceCollection &AutoTuner::getEvidenceCollection() const { return _evidenceCollection; }

bool AutoTuner::canMeasureEnergy() const { return _energyMeasurementPossible; }

void AutoTuner::setRebuildFrequency(double rebuildFrequency) { _rebuildFrequency = rebuildFrequency; }

void AutoTuner::checkEarlyStoppingCondition() {
  if (_samplesNotRebuildingNeighborLists.empty()) {
    // Wait until a sample without the expensive rebuilding has occurred to compare the configurations based on
    // estimated runtime for the non-tuning phase. This should generally happen in the second sample.
    return;
  }

  if (_evidenceCollection.empty()) {
    // Since there is no prior evidence, we must fully evaluate the first configuration.
    return;
  }

  long preliminaryEstimate = estimateRuntimeFromSamples();

  auto [_, bestEvidence] = _evidenceCollection.getLatestOptimalConfiguration();

  double slowdownFactor = static_cast<double>(preliminaryEstimate) / static_cast<double>(bestEvidence.value);

  if (slowdownFactor > _earlyStoppingFactor) {
    AutoPasLog(DEBUG,
               "Configuration is {} times slower than the current fastest traversal time. This is higher than the "
               "earlyStoppingFactor factor of {}. Further samples of this configuration will be skipped.",
               slowdownFactor, _earlyStoppingFactor);
    _earlyStoppingOfResampling = true;
  }
}

void AutoTuner::setContainerConstraint(ContainerOption container) { _containerConstraint = container; }

std::set<ContainerOption> AutoTuner::getSearchSpaceContainers() const {
  std::set<ContainerOption> containers;
  for (const auto &conf : _searchSpace) {
    containers.insert(conf.container);
  }
  return containers;
}

void AutoTuner::incrementTuningPhase() {
  ++_tuningPhase;
  if (_isTuning) {
    AutoPasLog(WARN, "Warning: Tuning needs more iterations than the specified tuning interval of {}!",
               _tuningInterval);
  }
  // New tuning is about to begin
  _requiresRebuild = true;
  // Incrementing the tuning phase means we are about to start tuning.
  _isTuning = true;
}
}  // namespace autopas
