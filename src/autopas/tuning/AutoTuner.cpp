/**
 * @file AutoTuner.cpp
 * @author F. Gratl
 * @date 21.06.23
 */

#include "AutoTuner.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <vector>

#include "autopas/tuning/selectors/OptimumSelector.h"
#include "autopas/tuning/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyLoggerWrapper.h"
#include "autopas/tuning/utils/Smoothing.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {
AutoTuner::AutoTuner(std::vector<std::unique_ptr<TuningStrategyInterface>> &tuningStrategies,
                     const std::set<Configuration> &searchSpace, const AutoTunerInfo &info)
    : _selectorStrategy(info.selectorStrategy),
      _tuningStrategies(info.useTuningStrategyLoggerProxy ? wrapTuningStrategies(tuningStrategies, info.outputSuffix)
                                                          : std::move(tuningStrategies)),
      _iteration(0),
      _tuningInterval(info.tuningInterval),
      _iterationsSinceTuning(info.tuningInterval),  // init to max so that tuning happens in first iteration
      _tuningMetric(info.tuningMetric),
      _energyMeasurementPossible(initEnergy()),
      _rebuildFrequency(info.rebuildFrequency),
      _maxSamples(info.maxSamples),
      _mpiTuningMaxDifferenceForBucket(info.MPITuningMaxDifferenceForBucket),
      _mpiTuningWeightForMaxDensity(info.MPITuningWeightForMaxDensity),
      _needsHomogeneityAndMaxDensity(std::transform_reduce(
          tuningStrategies.begin(), tuningStrategies.end(), false, std::logical_or(),
          [](auto &tuningStrat) { return tuningStrat->smoothedHomogeneityAndMaxDensityNeeded(); })),
      _needsLiveInfo(std::transform_reduce(tuningStrategies.begin(), tuningStrategies.end(), false, std::logical_or(),
                                           [](auto &tuningStrat) { return tuningStrat->needsLiveInfo(); })),
      _samplesNotRebuildingNeighborLists(info.maxSamples),
      _searchSpace(searchSpace),
      _configQueue(searchSpace.begin(), searchSpace.end()),
      _tuningResultLogger(info.outputSuffix),
      _tuningDataLogger(info.maxSamples, info.outputSuffix) {
  _homogeneitiesOfLastTenIterations.reserve(10);
  _maxDensitiesOfLastTenIterations.reserve(10);
  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
  }
  AutoPasLog(DEBUG, "Points in search space: {}", _searchSpace.size());
}

std::vector<std::unique_ptr<TuningStrategyInterface>> AutoTuner::wrapTuningStrategies(
    std::vector<std::unique_ptr<TuningStrategyInterface>> &tuningStrategies, const std::string &outputSuffix) {
  std::vector<std::unique_ptr<TuningStrategyInterface>> wrappedStrategies{};
  wrappedStrategies.reserve(tuningStrategies.size());
  std::transform(tuningStrategies.begin(), tuningStrategies.end(), std::back_inserter(wrappedStrategies),
                 [&](auto &tuningStrategy) {
                   return std::make_unique<TuningStrategyLoggerWrapper>(std::move(tuningStrategy), outputSuffix);
                 });
  return wrappedStrategies;
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

void AutoTuner::logIteration(const Configuration &conf, bool tuningIteration, long tuningTime) {
  // when a tuning result is found log it
  if (not tuningIteration) {
    AutoPasLog(DEBUG, "Selected Configuration {}", getCurrentConfig().toString());
    // TODO: When AutoTuner manages the search space also log smoothed time of the optimum.
    _tuningResultLogger.logTuningResult(conf, _iteration, tuningTime);
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

  // Determine where in a tuning phase we are
  if (_iterationsSinceTuning == _tuningInterval) {
    // CASE: Start of a tuning phase
    // in the first iteration of a tuning phase we reset all strategies
    // and refill the queue with the complete search space.
    std::copy(_searchSpace.begin(), _searchSpace.end(), std::back_inserter(_configQueue));
    // then let the strategies filter and sort it
    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
      tuningStrategy->reset(_iteration, _tuningPhase, _configQueue, _evidenceCollection);
    });
  } else if (_configQueue.empty()) {
    // CASE: End of a tuning phase
    // if the queue is empty we are done tuning.
    _iterationsSinceTuning = 0;
    const auto [optConf, optEvidence] = _evidenceCollection.getLatestOptimalConfiguration();
    _configQueue.push_back(optConf);
    stillTuning = false;
  } else {
    // CASE: somewhere in a tuning phase
    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
      return tuningStrategy->optimizeSuggestions(_configQueue, _evidenceCollection);
    });
    // samples are no longer needed.
    _samplesNotRebuildingNeighborLists.clear();
    _samplesRebuildingNeighborLists.clear();
  }
  tuningTimer.stop();

  return stillTuning;
}

const Configuration &AutoTuner::getCurrentConfig() const { return _configQueue.back(); }

std::tuple<Configuration, bool> AutoTuner::getNextConfig() {
  // If we are not (yet) tuning or there is nothing to tune return immediately.
  if (_iterationsSinceTuning < _tuningInterval or searchSpaceIsTrivial()) {
    return {getCurrentConfig(), false};
  } else if (getCurrentNumSamples() < _maxSamples) {
    // If we are still collecting samples from one config return immediately.
    return {getCurrentConfig(), true};
  } else {
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
  _configQueue.erase(
      std::remove_if(_configQueue.begin(), _configQueue.end(), [&](auto &conf) { return conf == rejectedConfig; }),
      _configQueue.end());

  if (indefinitely) {
    // delete rejected config from the search space and notify tuning strategies.
    _searchSpace.erase(rejectedConfig);
    std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(),
                  [&](auto &tuningStrategy) { tuningStrategy->rejectConfigurationIndefinitely(rejectedConfig); });
  }

  // let all configurations apply their optimizations in the order they are defined.
  // If any is still tuning consider the tuning phase still ongoing.
  std::for_each(_tuningStrategies.begin(), _tuningStrategies.end(), [&](auto &tuningStrategy) {
    return tuningStrategy->optimizeSuggestions(_configQueue, _evidenceCollection);
  });
  const auto stillTuning = _configQueue.empty();
  return {getCurrentConfig(), stillTuning};
}

void AutoTuner::addMeasurement(long sample, bool neighborListRebuilt) {
  const auto &currentConfig = _configQueue.back();
  if (getCurrentNumSamples() < _maxSamples) {
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
          std::min(reducedValue, smoothing::smoothLastPoint(_evidenceCollection.getEvidence(currentConfig), 5));

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
          }(),
          [&]() {
            std::ostringstream ss;
            // print config
            ss << currentConfig << " : ";
            // print all timings
            ss << utils::ArrayUtils::to_string(_samplesRebuildingNeighborLists, " ",
                                               {"With rebuilding neighbor lists [ ", " ]"});
            ss << utils::ArrayUtils::to_string(_samplesNotRebuildingNeighborLists, " ",
                                               {"Without rebuilding neighbor lists [ ", " ]"});
            ss << " Smoothed value: " << smoothedValue;
            return ss.str();
          }());

      _tuningDataLogger.logTuningData(currentConfig, _samplesRebuildingNeighborLists,
                                      _samplesNotRebuildingNeighborLists, _iteration, reducedValue, smoothedValue);

      // We finished processing this config so remove it from the queue
      _configQueue.pop_back();
    }
  }
}

void AutoTuner::bumpIterationCounters() {
  ++_iteration;
  ++_iterationsSinceTuning;
  // this will NOT catch the first tuning phase because _iterationsSinceTuning is initialized to _tuningInterval.
  // Hence, _tuningPhase is initialized as 1.
  if (_iterationsSinceTuning == _tuningInterval) {
    ++_tuningPhase;
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
  // Check if energy measurement is possible,
  try {
    _raplMeter.init();
    _raplMeter.reset();
    _raplMeter.sample();
  } catch (const utils::ExceptionHandler::AutoPasException &e) {
    if (_tuningMetric == TuningMetricOption::energy) {
      throw e;
    } else {
      AutoPasLog(WARN, "Energy Measurement not possible:\n\t{}", e.what());
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
        throw e;
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
        throw e;
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

  const auto numIterationsNotBuilding =
      std::max(0, static_cast<int>(_rebuildFrequency) - static_cast<int>(_samplesRebuildingNeighborLists.size()));
  const auto numIterationsBuilding = _rebuildFrequency - numIterationsNotBuilding;

  // calculate weighted estimate for one iteration
  return (numIterationsBuilding * reducedValueBuilding + numIterationsNotBuilding * reducedValueNotBuilding) /
         _rebuildFrequency;
}

std::tuple<bool, bool> AutoTuner::prepareIteration() {
  // Flag if this is the first iteration in a new tuning phase
  const bool startOfTuningPhase = _iterationsSinceTuning == _tuningInterval;

  // first tuning iteration -> reset everything
  if (startOfTuningPhase) {
    // reset the configQueue
    _configQueue.clear();
    _configQueue.reserve(_searchSpace.size());
    std::copy(_searchSpace.begin(), _searchSpace.end(), std::back_inserter(_configQueue));

    // call the appropriate versions of reset
    for (const auto &tuningStrat : _tuningStrategies) {
      if (auto *mpiStrategy = dynamic_cast<MPIParallelizedStrategy *>(tuningStrat.get())) {
        const std::pair<double, double> smoothedHomogeneityAndMaxDensity{
            autopas::OptimumSelector::medianValue(_homogeneitiesOfLastTenIterations),
            autopas::OptimumSelector::medianValue(_maxDensitiesOfLastTenIterations)};
        mpiStrategy->reset(_iteration, smoothedHomogeneityAndMaxDensity, _mpiTuningMaxDifferenceForBucket,
                           _mpiTuningWeightForMaxDensity);
      } else {
        tuningStrat->reset(_iteration, _tuningPhase, _configQueue, _evidenceCollection);
      }
    }

    // Homogeneity was calculated directly before the tuning phase so reset it now.
    if (_needsHomogeneityAndMaxDensity) {
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
  const bool needsLiveInfoNow = startOfTuningPhase and _needsLiveInfo;

  // calc homogeneity if needed, and we are within 10 iterations of the next tuning phase
  const size_t numIterationsForHomogeneity = 10;
  const bool needsHomogeneityAndMaxDensityNow =
      _needsHomogeneityAndMaxDensity and _iterationsSinceTuning > _tuningInterval - numIterationsForHomogeneity and
      _iterationsSinceTuning <= _tuningInterval;

  return {needsLiveInfoNow, needsHomogeneityAndMaxDensityNow};
}
}  // namespace autopas