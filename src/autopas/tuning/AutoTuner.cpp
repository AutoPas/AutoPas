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
      _tuningMetric(autoTunerInfo.tuningMetric),
      _useLOESSSmoothening(autoTunerInfo.useLOESSSmoothening),
      _energyMeasurementPossible(initEnergy()),
      _rebuildFrequency(rebuildFrequency),
      _maxSamples(autoTunerInfo.maxSamples),
      _earlyStoppingFactor(autoTunerInfo.earlyStoppingFactor),
      _needsDomainSimilarityStatistics(std::ranges::any_of(
          _tuningStrategies, [](auto &tuningStrat) { return tuningStrat->needsDomainSimilarityStatistics(); })),
      _needsLiveInfo(
          std::ranges::any_of(_tuningStrategies, [](auto &tuningStrat) { return tuningStrat->needsLiveInfo(); })),
      _samplesTraverseInteractions(autoTunerInfo.maxSamples),
      _searchSpace(searchSpace),
      _configQueue(searchSpace.begin(), searchSpace.end()),
      _tuningResultLogger(outputSuffix, autoTunerInfo.tuningMetric),
      _tuningDataLogger(autoTunerInfo.maxSamples, rebuildFrequency, outputSuffix),
      _energySensor(utils::EnergySensor(autoTunerInfo.energySensor)) {
  _samplesRebuildingNeighborLists.reserve(autoTunerInfo.maxSamples);
  _pdBinDensityStdDevOfLastTenIterations.reserve(10);
  _pdBinMaxDensityOfLastTenIterations.reserve(10);
  if (_searchSpace.empty()) {
    utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
  }
  AutoPasLog(DEBUG, "Points in search space: {}", _searchSpace.size());
  AutoPasLog(DEBUG, "AutoTuner constructed with LOESS Smoothening {}.", _useLOESSSmoothening ? "enabled" : "disabled");
}

AutoTuner &AutoTuner::operator=(AutoTuner &&other) noexcept {
  _tuningStrategies = std::move(other._tuningStrategies);
  return *this;
}

AutoTuner::~AutoTuner() = default;

void AutoTuner::addDomainSimilarityStatistics(double pdBinDensityStdDev, double pdBinMaxDensity) {
  _pdBinDensityStdDevOfLastTenIterations.push_back(pdBinDensityStdDev);
  _pdBinMaxDensityOfLastTenIterations.push_back(pdBinMaxDensity);
}

void AutoTuner::logTuningResult(const long tuningTime, const size_t currentIteration) const {
  // only log if we are at the end of a tuning phase
  if (_endOfTuningPhase) {
    // This string is part of several older scripts, hence it is not recommended to change it.
    const auto [conf, optimalEvidence] = _evidenceCollection.getLatestOptimalConfiguration();
    AutoPasLog(DEBUG, "Selected Configuration {}", conf.toString());
    _tuningResultLogger.logTuningResult(conf, currentIteration, tuningTime, optimalEvidence.reducedValue);
  }
}

bool AutoTuner::searchSpaceIsTrivial() const { return _searchSpace.size() == 1; }

bool AutoTuner::searchSpaceIsEmpty() const { return _searchSpace.empty(); }

void AutoTuner::forceRetune() {
  if (inTuningPhase()) {
    AutoPasLog(WARN, "Warning: Currently running tuning phase is aborted a new one is started!");
  }
  _samplesTraverseInteractions.resize(_maxSamples);
  _forceRetune = true;
}

bool AutoTuner::tuneConfiguration(const size_t currentIteration, const size_t tuningPhase,
                                  const bool isStartOfTuningPhase) {
  if (isStartOfTuningPhase or _forceRetune) {
    _isTuning = true;
  }
  // If we are not (yet) tuning or there is nothing to tune return immediately.
  if (not inTuningPhase()) {
    _isTuning = false;
    _forceRetune = false;
    return false;
  }
  if (getCurrentNumSamples() < _maxSamples and not _earlyStoppingOfResampling) {
    return true;
  }

  // We finished collection samples for this config so remove it from the queue
  if (not _configQueue.empty()) {
    _configQueue.pop_back();
  }

  // We plan to test a new config so clear all samples.
  _samplesTraverseInteractions.clear();
  _samplesRebuildingNeighborLists.clear();
  _earlyStoppingOfResampling = false;

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
  if (isStartOfTuningPhase or _forceRetune) {
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
    for (const auto &tuningStrategy : _tuningStrategies) {
      const auto configQueueBackup = _configQueue;
      const auto intentionalWipe =
          tuningStrategy->reset(currentIteration, tuningPhase, _configQueue, _evidenceCollection);
      AutoPasLog(DEBUG, "ConfigQueue after applying {}::reset(): (Size={}) {}",
                 tuningStrategy->getOptionType().to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
      if (not intentionalWipe) {
        restoreConfigQueueIfEmpty(configQueueBackup, tuningStrategy->getOptionType());
      }
    }
  } else {
    // CASE: somewhere in a tuning phase
    _isTuning = true;

    AutoPasLog(DEBUG, "ConfigQueue at tuneConfiguration before optimizeSuggestions: (Size={}) {}", _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
    for (const auto &tuningStrategy : _tuningStrategies) {
      const auto configQueueBackup = _configQueue;
      const auto intentionalWipe = tuningStrategy->optimizeSuggestions(_configQueue, _evidenceCollection);
      AutoPasLog(DEBUG, "ConfigQueue after applying {}::optimizeSuggestions(): (Size={}) {}",
                 tuningStrategy->getOptionType().to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
      if (not intentionalWipe) {
        restoreConfigQueueIfEmpty(configQueueBackup, tuningStrategy->getOptionType());
      }
    }
  }

  if (_configQueue.empty()) {
    // If the queue is empty we are done tuning.
    handleEndOfTuningPhase(tuningPhase);
  }
  return _isTuning;
}

void AutoTuner::handleEndOfTuningPhase(const size_t tuningPhase) {
  _endOfTuningPhase = true;
  _isTuning = false;
  selectBestConfiguration(tuningPhase);
}

void AutoTuner::selectBestConfiguration(const size_t tuningPhase) {
  // If there are still queued up configurations, clear them.
  _configQueue.clear();
  // Find and push_back the optimal configuration for the current container.
  const auto [optConf, optEvidence] =
      _evidenceCollection.getOptimalConfiguration(tuningPhase, EvidenceCollection::EvidenceMode::REDUCED, std::nullopt);
  _configQueue.push_back(optConf);

  _samplesRebuildingNeighborLists.resize(_maxSamples);
  _samplesTraverseInteractions.resize(_maxSamples);
}

void AutoTuner::forceOptimalConfiguration(const Configuration &optimalConfig) {
  _configQueue.clear();
  _configQueue.push_back(optimalConfig);

  // Fill up sample buffers to indicate we are not collecting samples anymore
  _samplesRebuildingNeighborLists.resize(_maxSamples);
  _samplesTraverseInteractions.resize(_maxSamples);

  _endOfTuningPhase = true;
  _isTuning = false;
}

const Configuration &AutoTuner::getCurrentConfig() const {
  if (_configQueue.empty()) {
    utils::ExceptionHandler::exception(
        "AutoTuner::getCurrentConfig(): Cannot get the current Configuration as the config queue is empty");
  }
  return _configQueue.back();
}

Configuration AutoTuner::rejectConfig(const Configuration &rejectedConfig, bool indefinitely,
                                      const size_t tuningPhase) {
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
  std::ranges::for_each(_tuningStrategies, [&](auto &tuningStrategy) {
    tuningStrategy->rejectConfiguration(rejectedConfig, indefinitely);
    AutoPasLog(DEBUG, "ConfigQueue after applying {}::rejectConfiguration(): (Size={}) {}",
               tuningStrategy->getOptionType().to_string(), _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
  });

  if (_configQueue.empty()) {
    handleEndOfTuningPhase(tuningPhase);
  }

  return getCurrentConfig();
}

void AutoTuner::addMeasurement(long sampleRebuild, long sampleTraverseParticles, const bool neighborListRebuilt,
                               const size_t iteration, const size_t tuningPhase) {
  const auto &currentConfig = _configQueue.back();
  AutoPasLog(TRACE, "Adding sampleRebuild and sampleNonRebuild {}, {} to configuration {}.", sampleRebuild,
             sampleTraverseParticles, currentConfig.toShortString());
  if (neighborListRebuilt) {
    // We add samples to _samplesRebuildingNeighborLists only for iterations where a neighbor list rebuild took place.
    // We do this to avoid essentially "zero" samples from the iterations without a neighbor list rebuild.
    // This is done so that we can apply different strategies (mean, median and min) to find the optimum value.
    _samplesRebuildingNeighborLists.push_back(sampleRebuild);
  }
  _samplesTraverseInteractions.push_back(sampleTraverseParticles);

  checkEarlyStoppingCondition();

  // if this was the last sample for this configuration:
  //  - calculate the evidence from the collected samples
  //  - log what was collected
  //  - remove the configuration from the queue
  if (getCurrentNumSamples() >= _maxSamples or _earlyStoppingOfResampling) {
    const long reducedValue = estimateRuntimeFromSamples();
    _evidenceCollection.addEvidence(currentConfig,
                                    {iteration, tuningPhase, reducedValue, sampleRebuild, sampleTraverseParticles});

    // If LOESS-based smoothening is enabled, use it to smooth evidence to remove high outliers. If smoothing results in
    // a higher value or if LOESS-based smoothening is disabled, use the original value.
    const auto smoothedValue =
        _useLOESSSmoothening
            ? std::min(reducedValue, smoothing::smoothLastPoint(*_evidenceCollection.getEvidence(currentConfig), 5))
            : reducedValue;

    // replace collected evidence with smoothed value to improve next smoothing
    _evidenceCollection.modifyLastEvidence(currentConfig).reducedValue = smoothedValue;

    for (const auto &tuningStrategy : _tuningStrategies) {
      tuningStrategy->addEvidence(getCurrentConfig(), _evidenceCollection.modifyLastEvidence(currentConfig));
    }

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
          ss << utils::ArrayUtils::to_string(_samplesTraverseInteractions, " ",
                                             {"Without rebuilding neighbor lists [ ", " ] "});
          ss << "Smoothed value: " << smoothedValue;
          return ss.str();
        }());

    auto samplesRebuildingNeighborLists = _samplesRebuildingNeighborLists;
    auto samplesTraverseInteractions = _samplesTraverseInteractions;

    // Pad the samples vectors  to the length _maxSamples as expected by the tuning data logger
    samplesRebuildingNeighborLists.resize(_maxSamples, -1);
    samplesTraverseInteractions.resize(_maxSamples, -1);

    _tuningDataLogger.logTuningData(currentConfig, samplesRebuildingNeighborLists, samplesTraverseInteractions,
                                    iteration, reducedValue, smoothedValue, _rebuildFrequency);
  }
}

bool AutoTuner::willRebuildNeighborLists() const {
  // This function should return true if the next call to tuneConfiguration() will lead to a new configuration being
  // selected. A new configuration will be selected if:
  //  The current configuration has finished collecting samples (either _maxSamples reached or early stopping).
  //  We prepare for a forced retuning.
  return _forceRetune or (_isTuning and (_earlyStoppingOfResampling or (getCurrentNumSamples() >= _maxSamples)));
}

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
  // We use the non rebuilding samples because these are added for every sample
  return _samplesTraverseInteractions.size();
}

long AutoTuner::estimateRuntimeFromSamples() const {
  // reduce samples for rebuild and non-rebuild iterations with the given selector strategy
  const auto reducedValueBuilding =
      autopas::OptimumSelector::optimumValue(_samplesRebuildingNeighborLists, _selectorStrategy);
  const auto reducedValueNotBuilding =
      autopas::OptimumSelector::optimumValue(_samplesTraverseInteractions, _selectorStrategy);

  // Calculate weighted average as if there was exactly one sample for each iteration in the rebuild interval.
  return reducedValueBuilding / _rebuildFrequency + reducedValueNotBuilding;
}

bool AutoTuner::needsLiveInfo() const { return (_needsLiveInfo or _needsDomainSimilarityStatistics); }

const std::vector<Configuration> &AutoTuner::getConfigQueue() const { return _configQueue; }

const std::vector<std::unique_ptr<TuningStrategyInterface>> &AutoTuner::getTuningStrategies() const {
  return _tuningStrategies;
}

void AutoTuner::receiveLiveInfo(const LiveInfo &liveInfo, const bool isStartOfTuningPhase) {
  // Handle Live Info processing before the tuning phase
  if (_needsDomainSimilarityStatistics) {
    const auto particleDependentBinDensityStdDev = liveInfo.get<double>("particleDependentBinDensityStdDev");
    const auto particleDependentBinMaxDensity = liveInfo.get<double>("particleDependentBinMaxDensity");
    addDomainSimilarityStatistics(particleDependentBinDensityStdDev, particleDependentBinMaxDensity);
  }
  // Handling at start of tuning phase
  if ((_needsLiveInfo or _needsDomainSimilarityStatistics) and isStartOfTuningPhase) {
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

    // Clear the vectors so the next tuning phase starts with fresh data
    _pdBinDensityStdDevOfLastTenIterations.clear();
    _pdBinMaxDensityOfLastTenIterations.clear();
  }
  if (_needsDomainSimilarityStatistics and _forceRetune and isStartOfTuningPhase) {
    // If we have sent statistics because of a forced retune, throw a warning
    AutoPasLog(WARN,
               "Due to a forced retune, domain similarity statistics, used by at least one tuning strategy, may"
               "not be correct.");
  }
}

const TuningMetricOption &AutoTuner::getTuningMetric() const { return _tuningMetric; }

bool AutoTuner::inTuningPhase() const {
  // If _iteration % _tuningInterval == 0 we are in the first tuning iteration but tuneConfiguration has not
  // been called yet.
  return (_isTuning or _forceRetune) and not searchSpaceIsTrivial();
}

// bool AutoTuner::inFirstTuningIteration() const { return (_iteration % _tuningInterval == 0); }

bool AutoTuner::inLastTuningIteration() const { return _endOfTuningPhase; }

bool AutoTuner::inFirstConfigurationLastSample() const {
  // todo: first config should be determined different in case a strategy filters out configs
  return (_configQueue.size() == _searchSpace.size()) and (getCurrentNumSamples() == _maxSamples - 1);
}

const EvidenceCollection &AutoTuner::getEvidenceCollection() const { return _evidenceCollection; }

bool AutoTuner::canMeasureEnergy() const { return _energyMeasurementPossible; }

void AutoTuner::setRebuildFrequency(double rebuildFrequency) { _rebuildFrequency = rebuildFrequency; }

void AutoTuner::checkEarlyStoppingCondition() {
  if (_evidenceCollection.empty()) {
    // Since there is no prior evidence, we must fully evaluate the first configuration.
    return;
  }

  long preliminaryEstimate = estimateRuntimeFromSamples();

  auto [_, bestEvidence] = _evidenceCollection.getLatestOptimalConfiguration();

  double slowdownFactor = static_cast<double>(preliminaryEstimate) / static_cast<double>(bestEvidence.reducedValue);

  if (slowdownFactor > _earlyStoppingFactor) {
    AutoPasLog(DEBUG,
               "Configuration is {} times slower than the current fastest traversal time. This is higher than the "
               "earlyStoppingFactor factor of {}. Further samples of this configuration will be skipped.",
               slowdownFactor, _earlyStoppingFactor);
    _earlyStoppingOfResampling = true;
  }
}

std::set<ContainerOption> AutoTuner::getSearchSpaceContainers() const {
  std::set<ContainerOption> containers;
  for (const auto &conf : _searchSpace) {
    containers.insert(conf.container);
  }
  return containers;
}
}  // namespace autopas
