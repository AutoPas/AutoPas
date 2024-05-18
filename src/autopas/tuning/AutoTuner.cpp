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
      _iterationsSinceTuning(autoTunerInfo.tuningInterval),  // init to max so that tuning happens in first iteration
      _stillTuning(searchSpace.size() > 1),
      _tuningMetric(autoTunerInfo.tuningMetric),
      _energyMeasurementPossible(initEnergy()),
      _rebuildFrequency(rebuildFrequency),
      _maxSamples(autoTunerInfo.maxSamples),
      _needsHomogeneityAndMaxDensity(std::transform_reduce(
          _tuningStrategies.begin(), _tuningStrategies.end(), false, std::logical_or(),
          [](auto &tuningStrat) { return tuningStrat->needsSmoothedHomogeneityAndMaxDensity(); })),
      _needsLiveInfo(std::transform_reduce(_tuningStrategies.begin(), _tuningStrategies.end(), false, std::logical_or(),
                                           [](auto &tuningStrat) { return tuningStrat->needsLiveInfo(); })),
      _samplesNotRebuildingNeighborLists(autoTunerInfo.maxSamples),
      _searchSpace(searchSpace),
      _configQueue(searchSpace.begin(), searchSpace.end()),
      _tuningResultLogger(outputSuffix),
      _tuningDataLogger(autoTunerInfo.maxSamples, outputSuffix) {
  _samplesRebuildingNeighborLists.reserve(autoTunerInfo.maxSamples);
  _homogeneitiesOfLastTenIterations.reserve(10);
  _maxDensitiesOfLastTenIterations.reserve(10);
  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
  }
  AutoPasLog(DEBUG, "Points in search space: {}", _searchSpace.size());
}

AutoTuner &AutoTuner::operator=(AutoTuner &&other) noexcept {
  _tuningStrategies = std::move(other._tuningStrategies);
  return *this;
}

void AutoTuner::addHomogeneityAndMaxDensity(double homogeneity, double maxDensity, long time) {
  _homogeneitiesOfLastTenIterations.push_back(homogeneity);
  _maxDensitiesOfLastTenIterations.push_back(maxDensity);
  _timerCalculateHomogeneity.addTime(time);
}

void AutoTuner::logTuningResult(bool tuningIteration, long tuningTime) {
  // only log if we are at the end of a tuning phase
  if (not tuningIteration and _iterationsSinceTuning == 0 and not _evidenceCollection.empty()) {
    // This string is part of several older scripts, hence it is not recommended to change it.
    const auto [conf, optimalEvidence] = _evidenceCollection.getLatestOptimalConfiguration();
    AutoPasLog(DEBUG, "Selected Configuration {}", conf.toString());
    _tuningResultLogger.logTuningResult(conf, _iteration, tuningTime, optimalEvidence.value);
  }
}

bool AutoTuner::searchSpaceIsTrivial() const { return _searchSpace.size() == 1; }

bool AutoTuner::searchSpaceIsEmpty() const { return _searchSpace.empty(); }

void AutoTuner::forceRetune() {
  _iterationsSinceTuning = _tuningInterval;
  _samplesNotRebuildingNeighborLists.resize(_maxSamples);
}

bool AutoTuner::tuneConfiguration() {
  bool stillTuning = true;

  utils::Timer tuningTimer;
  tuningTimer.start();

  // We finished collection samples for this config so remove it from the queue
  _configQueue.pop_back();

  // We plan to test a new config so clear all samples.
  _samplesNotRebuildingNeighborLists.clear();
  _samplesRebuildingNeighborLists.clear();

  // Helper function to reset the ConfigQueue if something wipes it.
  auto restoreConfigQueueIfEmpty = [&](const auto &configQueueBackup, const TuningStrategyOption &stratOpt) {
    if (_configQueue.empty()) {
      _configQueue = configQueueBackup;
    }
    AutoPasLog(WARN, "ConfigQueue wipe by {} detected! Resetting to previous state: (Size={}) {}", stratOpt.to_string(),
               _configQueue.size(), utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"}, [](const auto &conf) {
                 return conf.toShortString(false);
               }));
  };

  // Determine where in a tuning phase we are
  if (_iterationsSinceTuning == _tuningInterval) {
    // CASE: Start of a tuning phase
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
      tuningStrategy->reset(_iteration, _tuningPhase, _configQueue, _evidenceCollection);
      AutoPasLog(DEBUG, "ConfigQueue after applying {}::reset(): (Size={}) {}",
                 tuningStrategy->getOptionType().to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
      restoreConfigQueueIfEmpty(configQueueBackup, tuningStrategy->getOptionType());
    });
  } else {
    // CASE: somewhere in a tuning phase
    AutoPasLog(DEBUG, "ConfigQueue at tuneConfiguration before optimizeSuggestions: (Size={}) {}", _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
      const auto configQueueBackup = _configQueue;
      tuningStrategy->optimizeSuggestions(_configQueue, _evidenceCollection);
      AutoPasLog(DEBUG, "ConfigQueue after applying {}::optimizeSuggestions(): (Size={}) {}",
                 tuningStrategy->getOptionType().to_string(), _configQueue.size(),
                 utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                              [](const auto &conf) { return conf.toShortString(false); }));
      restoreConfigQueueIfEmpty(configQueueBackup, tuningStrategy->getOptionType());
    });
  }

  // CASE: End of a tuning phase. This is not exclusive to the other cases!
  if (_configQueue.empty()) {
    const auto [optConf, optEvidence] = _evidenceCollection.getOptimalConfiguration(_tuningPhase);
    _configQueue.push_back(optConf);
    stillTuning = false;
    // fill up sample buffer to indicate we are not collecting samples anymore
    _samplesRebuildingNeighborLists.resize(_maxSamples);
  }
  tuningTimer.stop();
  _stillTuning = stillTuning;

  return stillTuning;
}

const Configuration &AutoTuner::getCurrentConfig() const { return _configQueue.back(); }

std::tuple<Configuration, bool> AutoTuner::getNextConfig() {
  // If we are not (yet) tuning or there is nothing to tune return immediately.
  if (not inTuningPhase()) {
    return {getCurrentConfig(), false};
  } else if (getCurrentNumSamples() < _maxSamples) {
    // If we are still collecting samples from one config return immediately.
    return {getCurrentConfig(), true};
  } else {
    // This case covers any iteration in a tuning phase where a new configuration is needed (even the start of a phase)
    // If we are at the start of a phase tuneConfiguration() will also refill the queue and call reset on all strategies
    const bool stillTuning = tuneConfiguration();
    return {getCurrentConfig(), stillTuning};
  }
}

std::tuple<Configuration, bool> AutoTuner::rejectConfig(const Configuration &rejectedConfig, bool indefinitely) {
  if (searchSpaceIsTrivial()) {
    utils::ExceptionHandler::exception("Rejected the only configuration in the search space!\n{}",
                                       rejectedConfig.toString());
  }

  // remove the config from the queue
  _configQueue.erase(std::remove_if(_configQueue.begin(), _configQueue.end(),
                                    [&](const auto &conf) { return conf == rejectedConfig; }),
                     _configQueue.end());

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

  // let all strategies optimize the queue in the order they are defined.
  // If any is still tuning consider the tuning phase still ongoing.
  std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
    tuningStrategy->optimizeSuggestions(_configQueue, _evidenceCollection);
    AutoPasLog(DEBUG, "ConfigQueue after applying {}::optimizeSuggestions(): (Size={}) {}",
               tuningStrategy->getOptionType().to_string(), _configQueue.size(),
               utils::ArrayUtils::to_string(_configQueue, ", ", {"[", "]"},
                                            [](const auto &conf) { return conf.toShortString(false); }));
  });
  const auto stillTuning = not _configQueue.empty();
  return {getCurrentConfig(), stillTuning};
}

void AutoTuner::addMeasurement(long sample, bool neighborListRebuilt) {
  const auto &currentConfig = _configQueue.back();
  // sanity check
  if (getCurrentNumSamples() >= _maxSamples) {
    utils::ExceptionHandler::exception(
        "AutoTuner::addMeasurement(): Trying to add a new measurement to the AutoTuner but there are already enough"
        "for this configuration!\n"
        "tuneConfiguration() should have been called before to process and flush samples.");
  }
  AutoPasLog(TRACE, "Adding sample {} to configuration {}.", sample, currentConfig.toShortString());
  if (neighborListRebuilt) {
    _samplesRebuildingNeighborLists.push_back(sample);
  } else {
    _samplesNotRebuildingNeighborLists.push_back(sample);
  }
  // if this was the last sample for this configuration:
  //  - calculate the evidence from the collected samples
  //  - log what was collected
  //  - remove the configuration from the queue
  if (getCurrentNumSamples() == _maxSamples) {
    const long reducedValue = estimateRuntimeFromSamples();
    _evidenceCollection.addEvidence(currentConfig, {_iteration, _tuningPhase, reducedValue});

    // smooth evidence to remove high outliers. If smoothing results in a higher value use the original value.
    const auto smoothedValue =
        std::min(reducedValue, smoothing::smoothLastPoint(*_evidenceCollection.getEvidence(currentConfig), 5));

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

    _tuningDataLogger.logTuningData(currentConfig, _samplesRebuildingNeighborLists, _samplesNotRebuildingNeighborLists,
                                    _iteration, reducedValue, smoothedValue);
  }
}

void AutoTuner::bumpIterationCounters(bool needToWait) {
  // reset counter after all autotuners finished tuning
  if (not(needToWait or inTuningPhase() or _iterationsSinceTuning < _tuningInterval)) {
    _iterationsSinceTuning = 0;
  }
  ++_iterationsSinceTuning;
  ++_iteration;

  // this will NOT catch the first tuning phase because _iterationsSinceTuning is initialized to _tuningInterval.
  // Hence, _tuningPhase is initialized as 1.
  if (_iterationsSinceTuning == _tuningInterval) {
    ++_tuningPhase;
    _stillTuning = true;
  }
}

bool AutoTuner::willRebuildNeighborLists() const {
  const bool inTuningPhase = _iterationsSinceTuning >= _tuningInterval;
  // How many iterations ago did the rhythm of rebuilds change?
  const auto iterationBaseline = inTuningPhase ? (_iterationsSinceTuning - _tuningInterval) : _iterationsSinceTuning;
  // What is the rebuild rhythm?
  const auto iterationsPerRebuild = inTuningPhase ? _maxSamples : _rebuildFrequency;
  return (iterationBaseline % iterationsPerRebuild) == 0;
}

bool AutoTuner::initEnergy() {
  // Check if energy measurement is possible.
  std::string errMsg(_raplMeter.init());
  if (errMsg.empty()) {
    try {
      _raplMeter.reset();
      _raplMeter.sample();
    } catch (const utils::ExceptionHandler::AutoPasException &e) {
      if (_tuningMetric == TuningMetricOption::energy) {
        utils::ExceptionHandler::exception(e);
      } else {
        AutoPasLog(WARN, "Energy Measurement not possible:\n\t{}", e.what());
        return false;
      }
    }
  } else {
    if (_tuningMetric == TuningMetricOption::energy) {
      utils::ExceptionHandler::exception(errMsg);
    } else {
      AutoPasLog(WARN, "Energy Measurement not possible:\n\t{}", errMsg);
      return false;
    }
  }
  return true;
}

bool AutoTuner::resetEnergy() {
  if (_energyMeasurementPossible) {
    try {
      _raplMeter.reset();
    } catch (const utils::ExceptionHandler::AutoPasException &e) {
      /**
       * very unlikely to happen, as check was performed at initialisation of autotuner
       * but may occur if permissions are changed during runtime.
       */
      AutoPasLog(WARN, "Energy Measurement no longer possible:\n\t{}", e.what());
      _energyMeasurementPossible = false;
      if (_tuningMetric == TuningMetricOption::energy) {
        utils::ExceptionHandler::exception(e);
      }
    }
  }
  return _energyMeasurementPossible;
}

std::tuple<double, double, double, long> AutoTuner::sampleEnergy() {
  if (_energyMeasurementPossible) {
    try {
      _raplMeter.sample();
    } catch (const utils::ExceptionHandler::AutoPasException &e) {
      AutoPasLog(WARN, "Energy Measurement no longer possible:\n\t{}", e.what());
      _energyMeasurementPossible = false;
      if (_tuningMetric == TuningMetricOption::energy) {
        utils::ExceptionHandler::exception(e);
      }
    }
  }
  return {_raplMeter.get_psys_energy(), _raplMeter.get_pkg_energy(), _raplMeter.get_ram_energy(),
          _raplMeter.get_total_energy()};
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

bool AutoTuner::prepareIteration() {
  // Flag if this is the first iteration in a new tuning phase
  const bool startOfTuningPhase = _iterationsSinceTuning == _tuningInterval;

  // first tuning iteration -> reset everything
  if (startOfTuningPhase) {
    // If needed, calculate homogeneity and maxDensity, and reset buffers.
    const auto [homogeneity, maxDensity] = [&]() {
      if (_needsHomogeneityAndMaxDensity) {
        const auto retTuple = std::make_tuple(OptimumSelector::medianValue(_homogeneitiesOfLastTenIterations),
                                              OptimumSelector::medianValue(_maxDensitiesOfLastTenIterations));
        _homogeneitiesOfLastTenIterations.clear();
        _maxDensitiesOfLastTenIterations.clear();
        AutoPasLog(DEBUG, "Calculating homogeneities over 10 iterations took in total {} ns on rank {}.",
                   _timerCalculateHomogeneity.getTotalTime(), []() {
                     int rank{0};
                     AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
                     return rank;
                   });
        return retTuple;
      } else {
        return std::make_tuple(-1., -1.);
      }
    }();

    // pass homogeneity and maxDensity info if needed
    for (const auto &tuningStrat : _tuningStrategies) {
      tuningStrat->receiveSmoothedHomogeneityAndMaxDensity(homogeneity, maxDensity);
    }
  }

  // if necessary, we need to collect live info in the first tuning iteration
  const bool needsLiveInfoNow = startOfTuningPhase and _needsLiveInfo;

  return needsLiveInfoNow;
}

bool AutoTuner::needsHomogeneityAndMaxDensityBeforePrepare() const {
  // calc homogeneity if needed, and we are within 10 iterations of the next tuning phase
  constexpr size_t numIterationsForHomogeneity = 10;
  return _needsHomogeneityAndMaxDensity and _iterationsSinceTuning > _tuningInterval - numIterationsForHomogeneity and
         _iterationsSinceTuning <= _tuningInterval;
}

const std::vector<Configuration> &AutoTuner::getConfigQueue() const { return _configQueue; }

const std::vector<std::unique_ptr<TuningStrategyInterface>> &AutoTuner::getTuningStrategies() const {
  return _tuningStrategies;
}

void AutoTuner::receiveLiveInfo(const LiveInfo &liveInfo) {
  for (auto &tuningStrategy : _tuningStrategies) {
    tuningStrategy->receiveLiveInfo(liveInfo);
  }
}

const TuningMetricOption &AutoTuner::getTuningMetric() const { return _tuningMetric; }

bool AutoTuner::inTuningPhase() const { return _stillTuning and not searchSpaceIsTrivial(); }

const EvidenceCollection &AutoTuner::getEvidenceCollection() const { return _evidenceCollection; }
}  // namespace autopas
