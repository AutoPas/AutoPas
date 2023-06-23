/**
 * @file AutoTuner.cpp
 * @author F. Gratl
 * @date 21.06.23
 */

#include "AutoTuner.h"

#include "autopas/selectors/Smoothing.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyLoggerWrapper.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {
AutoTuner::AutoTuner(std::unique_ptr<TuningStrategyInterface> tuningStrategy, double MPITuningMaxDifferenceForBucket,
                     double MPITuningWeightForMaxDensity, SelectorStrategyOption selectorStrategy,
                     TuningMetricOption tuningMetric, unsigned int tuningInterval, unsigned int maxSamples,
                     unsigned int rebuildFrequency, const std::string &outputSuffix, bool useTuningStrategyLoggerProxy)
    : _selectorStrategy(selectorStrategy),
      _tuningStrategy(useTuningStrategyLoggerProxy
                          ? std::make_unique<TuningStrategyLoggerWrapper>(std::move(tuningStrategy), outputSuffix)
                          : std::move(tuningStrategy)),
      _iteration(0),
      _tuningInterval(tuningInterval),
      _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
      _tuningMetric(tuningMetric),
      _energyMeasurementPossible(initEnergy()),
      _rebuildFrequency(rebuildFrequency),
      _maxSamples(maxSamples),
      _mpiTuningMaxDifferenceForBucket(MPITuningMaxDifferenceForBucket),
      _mpiTuningWeightForMaxDensity(MPITuningWeightForMaxDensity),
      _samplesNotRebuildingNeighborLists(maxSamples),
      _tuningResultLogger(outputSuffix),
      _tuningDataLogger(maxSamples, outputSuffix) {
  _homogeneitiesOfLastTenIterations.reserve(10);
  _maxDensitiesOfLastTenIterations.reserve(10);
  if (_tuningStrategy->searchSpaceIsEmpty()) {
    autopas::utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
  }
  AutoPasLog(DEBUG, "Points in search space: {}", _searchSpace.size());
}

AutoTuner &AutoTuner::operator=(AutoTuner &&other) noexcept {
  _tuningStrategy = std::move(other._tuningStrategy);
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

/**
 * Access to the searchSpaceIsTrivial bool variable (true if search space size  is 1 or less).
 * @return Smart pointer to the searchSpaceIsTrivial variable.
 */

bool AutoTuner::searchSpaceIsTrivial() { return _tuningStrategy->searchSpaceIsTrivial(); }

void AutoTuner::forceRetune() {
  _iterationsSinceTuning = _tuningInterval;
  _samplesNotRebuildingNeighborLists.resize(_maxSamples);
}

bool AutoTuner::tune() {
  bool stillTuning = true;

  utils::Timer tuningTimer;
  tuningTimer.start();

  // in the first iteration of a tuning phase we go with the initial state of the strategy.
  // Hence, only tune if we are not in the first iteration.
  if (_iterationsSinceTuning != _tuningInterval) {
    stillTuning = _tuningStrategy->tune(false);
  }

  // If we reach this line and are still tuning we have a new config, hence, we need to clear the samples.
  if (stillTuning) {
    // samples are no longer needed.
    _samplesNotRebuildingNeighborLists.clear();
    _samplesRebuildingNeighborLists.clear();
  }
  tuningTimer.stop();

  return stillTuning;
}

const Configuration &AutoTuner::getCurrentConfig() const { return _tuningStrategy->getCurrentConfiguration(); }

std::tuple<Configuration, bool> AutoTuner::getNextConfig() {
  // If we are not (yet) tuning or there is nothing to tune return immediately.
  if (_iterationsSinceTuning < _tuningInterval or searchSpaceIsTrivial()) {
    return {getCurrentConfig(), false};
  } else if (getCurrentNumSamples() < _maxSamples) {
    // If we are still collecting samples from one config return immediately.
    return {getCurrentConfig(), true};
  } else {
    const bool stillTuning = tune();
    return {getCurrentConfig(), stillTuning};
  }
}

void AutoTuner::addMeasurement(long sample, bool neighborListRebuilt) {
  const auto &currentConfig = _tuningStrategy->getCurrentConfiguration();
  if (getCurrentNumSamples() < _maxSamples) {
    AutoPasLog(TRACE, "Adding sample.");
    if (neighborListRebuilt) {
      _samplesRebuildingNeighborLists.push_back(sample);
    } else {
      _samplesNotRebuildingNeighborLists.push_back(sample);
    }
    // if this was the last sample:
    if (getCurrentNumSamples() == _maxSamples) {
      auto &evidenceCurrentConfig = _searchSpace[currentConfig];

      const long reducedValue = estimateRuntimeFromSamples();

      evidenceCurrentConfig.emplace_back(_iteration, reducedValue);

      // smooth evidence to remove high outliers. If smoothing results in a higher value use the original value.
      const auto smoothedValue = std::min(reducedValue, smoothing::smoothLastPoint(evidenceCurrentConfig, 5));

      // replace collected evidence with smoothed value to improve next smoothing
      evidenceCurrentConfig.back().value = smoothedValue;

      _tuningStrategy->addEvidence(smoothedValue, _iteration);

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
    }
  }
}

void AutoTuner::bumpIterationCounters() {
  ++_iteration;
  ++_iterationsSinceTuning;
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

}  // namespace autopas