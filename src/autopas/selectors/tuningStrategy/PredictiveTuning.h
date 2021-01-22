/**
 * @file PredictiveTuning.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <limits>
#include <set>
#include <utility>

#include "SetSearchSpaceBasedTuningStrategy.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * On every tuning phase, this strategy makes a runtime prediction for every configuration and then
 * only tests those which are within a certain range of the best prediction. In the end, the configuration
 * that performed best during testing is selected.
 */
class PredictiveTuning : public SetSearchSpaceBasedTuningStrategy {
 public:
  /**
   * Constructor for the PredictiveTuning that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param relativeOptimum
   * @param maxTuningIterationsWithoutTest
   * @param relativeRangeForBlacklist
   * @param testsUntilFirstPrediction
   * @param extrapolationMethodOption
   */
  PredictiveTuning(const std::set<ContainerOption> &allowedContainerOptions,
                   const std::set<double> &allowedCellSizeFactors,
                   const std::set<TraversalOption> &allowedTraversalOptions,
                   const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                   const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                   const std::set<Newton3Option> &allowedNewton3Options, double relativeOptimum,
                   unsigned int maxTuningIterationsWithoutTest, double relativeRangeForBlacklist,
                   unsigned int testsUntilFirstPrediction, ExtrapolationMethodOption extrapolationMethodOption)
      : SetSearchSpaceBasedTuningStrategy(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                                          allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options),
        _notTestedYet(_searchSpace),
        _relativeOptimumRange(relativeOptimum),
        _maxTuningIterationsWithoutTest(maxTuningIterationsWithoutTest),
        _relativeBlacklistRange(relativeRangeForBlacklist),
        _extrapolationMethod(extrapolationMethodOption) {
    setTestsUntilFirstPrediction(testsUntilFirstPrediction);

    // sets traversalTimesStorage
    for (const auto &configuration : _searchSpace) {
      std::vector<std::pair<size_t, long>> vector;
      _traversalTimesStorage.emplace(configuration, vector);
    }
  }

  /**
   * Constructor for the PredictiveTuning that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit PredictiveTuning(std::set<Configuration> allowedConfigurations)
      : SetSearchSpaceBasedTuningStrategy(std::move(allowedConfigurations)), _notTestedYet(_searchSpace) {}

  inline void addEvidence(long time, size_t iteration) override {
    _traversalTimesStorage[*_currentConfig].emplace_back(iteration, time);
    _lastTest[*_currentConfig] = _tuningPhaseCounter;
  }

  inline long getEvidence(Configuration configuration) const override {
    // compute the average of times for this configuration
    auto times = _traversalTimesStorage.at(configuration);
    long result = 0;
    for (auto time : times) {
      result += time.second;
    }
    return result / times.size();
  }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void reset(size_t iteration) override {
    _configurationPredictions.clear();
    _optimalSearchSpace.clear();
    _tooLongNotTestedSearchSpace.clear();
    _validSearchSpace = _searchSpace;
    _validConfigurationFound = false;
    _iterationBeginTuningPhase = iteration;

    selectOptimalSearchSpace();
  }

  inline bool tune(bool currentInvalid = false) override;

 private:
  /**
   * Set the value of _evidenceFirstPrediction
   * @param testsUntilFirstPrediction
   */
  inline void setTestsUntilFirstPrediction(unsigned int testsUntilFirstPrediction);

  inline void selectOptimalConfiguration();
  /**
   * Selects the configurations that are going to be tested.
   */
  inline void selectOptimalSearchSpace();
  /**
   * Provides different extrapolation methods for the prediction of the traversal time.
   */
  inline void calculatePredictions();
  /**
   * Predicts the traversal time by placing a line through the last two traversal points and calculating the prediction
   * for the current time.
   */
  inline void linePrediction();
  /**
   * Predicts the traversal time by creating a function that places a line through the data points and calculating the
   * prediction for the current time.
   */
  inline void linearRegression();
  /**
   * Creates a polynomial function using the Lagrange interpolation and with this function the prediction is calculated.
   */
  inline void lagrangePolynomial();
  /**
   * Creates a polynomial function using Newton's method of finite differences and with this function the prediction is
   * calculated.
   */
  inline void newtonPolynomial();
  /**
   * Creates a new optimalSearchSpace if every configuration in the previous one was invalid.
   */
  inline void reselectOptimalSearchSpace();
  /**
   * Selects configurations that go on the testing black list.
   */
  inline void blacklistBadConfigurations();
  /**
   * Creates output of the configuration and the prediction.
   */
  inline void predictionOutput(Configuration configuration);

  /**
   * Stores the traversal times for each configuration.
   * @param Configuration
   * @param Vector with pairs of the iteration and the traversal time.
   */
  std::unordered_map<Configuration, std::vector<std::pair<size_t, long>>, ConfigHash> _traversalTimesStorage;
  /**
   * A Map that for each configuration stores the function for the prediction to reuse it if no new traversal time was
   * added in the last tuning  phase. The way that function is stored depends on the prediction method, hence it is a
   * vector. Line Prediction: Gradient and last evidence Linear Regression: Gradient and iteration Newton: Vector of
   * coefficients
   */
  std::unordered_map<Configuration, std::vector<double>, ConfigHash> _configurationPredictionFunction;

  /**
   * Contains the predicted time for each configuration.
   * @param Configuration
   * @param traversal prediction
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _configurationPredictions;
  /**
   * Contains the configuration that are predicted to be optimal and going to be tested.
   */
  std::set<Configuration> _optimalSearchSpace;
  /**
   * Contains the configuration that have not been tested for a period of time and are going to be tested.
   */
  std::set<Configuration> _tooLongNotTestedSearchSpace;
  /**
   * Contains the configurations that are not marked invalid in the current tuning phase.
   */
  std::set<Configuration> _validSearchSpace;
  /**
   * Contains the configurations that have not been tested yet.
   */
  std::set<Configuration> _notTestedYet;
  /**
   * Stores the the last tuning phase a configuration got tested.
   * @param Configuration
   * @param last tuning phase
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _lastTest;
  /**
   * Gets incremented after every completed tuning phase.
   * This is initialized with one, because in selectOptimalConfiguration() there is a comparison with _lastTest and in
   * the first tuning phase it would use every configuration, if it was tested or not, because the _lastTest would
   * return zero in both cases.
   */
  unsigned int _tuningPhaseCounter{1};
  /**
   * Stores the iteration at the beginning of a tuning phase.
   * Mainly used for the traversal time storage.
   */
  unsigned int _iterationBeginTuningPhase{0};
  /**
   * Indicates if a valid configuration was found in the _optimalSearchSpace.
   */
  bool _validConfigurationFound = false;
  /**
   * Factor of the range of the optimal configurations for the optimalSearchSpace.
   */
  const double _relativeOptimumRange{1.2};
  /**
   * After not being tested this number of tuningPhases a configuration is being emplaced in _optimalSearchSpace.
   */
  const unsigned int _maxTuningIterationsWithoutTest{5};
  /**
   * The relative cutoff for configurations to be blacklisted.
   */
  const double _relativeBlacklistRange{0};
  /**
   * Stores the extrapolation method that is going to be used for the traversal time predictions.
   */
  const ExtrapolationMethodOption _extrapolationMethod{ExtrapolationMethodOption::linearRegression};
  /**
   * Stores the number of tests that have to be made until the first prediction.
   * This number also determines how much evidence is used for the prediction and for a polynomial extrapolation method,
   * this is also the degree of the polynomial.
   */
  unsigned int _evidenceFirstPrediction{3};
};

void PredictiveTuning::setTestsUntilFirstPrediction(unsigned int testsUntilFirstPrediction) {
  switch (_extrapolationMethod) {
    case ExtrapolationMethodOption::linePrediction: {
      _evidenceFirstPrediction = 2;
      break;
    }
    default: {
      _evidenceFirstPrediction = testsUntilFirstPrediction;
    }
  }
}

void PredictiveTuning::selectOptimalSearchSpace() {
  // Special case: Search space is trivial or more evidence is needed.
  if (_searchSpace.size() == 1 or _tuningPhaseCounter < _evidenceFirstPrediction + 1) {
    if (_searchSpace.empty()) {
      autopas::utils::ExceptionHandler::exception("PredictiveTuning: No possible configuration prediction found!");
    }
    _currentConfig = _searchSpace.begin();
    return;
  }

  calculatePredictions();

  const auto optimum = getOptimum(_configurationPredictions);

  _optimalSearchSpace.emplace(optimum->first);

  // selects configurations that are near the optimal prediction or have not been tested in a certain number of
  // iterations
  for (const auto &configuration : _searchSpace) {
    // Adds configurations that have not been tested for _maxTuningIterationsWithoutTest or are within the
    // _relativeOptimumRange
    if (static_cast<double>(_configurationPredictions[configuration]) / optimum->second <= _relativeOptimumRange) {
      _optimalSearchSpace.emplace(configuration);
    } else if (_tuningPhaseCounter - _lastTest[configuration] > _maxTuningIterationsWithoutTest) {
      _tooLongNotTestedSearchSpace.emplace(configuration);
    }
  }

  // sanity check
  if (_optimalSearchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("PredictiveTuning: No possible configuration prediction found!");
  }

  _currentConfig = _optimalSearchSpace.begin();
}

void PredictiveTuning::calculatePredictions() {
  switch (_extrapolationMethod) {
    case ExtrapolationMethodOption::linePrediction: {
      linePrediction();
      break;
    }
    case ExtrapolationMethodOption::linearRegression: {
      linearRegression();
      break;
    }
    case ExtrapolationMethodOption::lagrange: {
      lagrangePolynomial();
      break;
    }
    case ExtrapolationMethodOption::newton: {
      newtonPolynomial();
      break;
    }
  }
}

void PredictiveTuning::linePrediction() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    // if configuration was not tested in last tuning phase reuse prediction function.
    // checks is _configurationPredictionFunctionParams has a linear function saved and only then it can reuse it
    if (_lastTest[configuration] != (_tuningPhaseCounter - 1) and functionParams.size() == 2) {
      const auto delta = _iterationBeginTuningPhase - _traversalTimesStorage[configuration].back().first;

      // gradient * delta + last point
      _configurationPredictions[configuration] = functionParams[0] * delta + functionParams[1];
      logPrediction(configuration);
    } else if (const auto &traversalValues = _traversalTimesStorage[configuration];
               traversalValues.size() >= _evidenceFirstPrediction) {
      const auto &[traversal1Iteration, traversal1Time] = traversalValues[traversalValues.size() - 1];
      const auto &[traversal2Iteration, traversal2Time] = traversalValues[traversalValues.size() - 2];

      const auto gradient = (traversal1Time - traversal2Time) / (traversal1Iteration - traversal2Iteration);
      const auto delta = _iterationBeginTuningPhase - traversal1Iteration;

      // time1 + (time1 - time2) / (iteration1 - iteration2) / tuningPhase - iteration1)
      _configurationPredictions[configuration] = traversal1Time + gradient * delta;
      functionParams.clear();
      functionParams.emplace_back(gradient);
      functionParams.emplace_back(traversal1Time);
      logPrediction(configuration);
    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = std::numeric_limits<long unsigned int>::max();
      _tooLongNotTestedSearchSpace.emplace(configuration);
      AutoPasLog(debug, "No traversal time prediction for {}", configuration.toString());
    }
  }
}

void PredictiveTuning::linearRegression() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    // if configuration was not tested in last tuning phase reuse prediction function.
    // checks is _configurationPredictionFunctionParams has a linear function saved and only then it can reuse it
    if ((_lastTest[configuration] != (_tuningPhaseCounter - 1)) and functionParams.size() == 2) {
      // if configuration was not tested in last tuning phase reuse prediction function.
      // gradient * iteration + y-intercept
      _configurationPredictions[configuration] = functionParams[0] * _iterationBeginTuningPhase + functionParams[1];
      logPrediction(configuration);
    } else if (const auto &traversalValues = _traversalTimesStorage[configuration];
               traversalValues.size() >= _evidenceFirstPrediction) {
      size_t iterationMultTime = 0, iterationSum = 0, iterationSquareSum = 0, timeSum = 0;

      for (auto i = traversalValues.size() - _evidenceFirstPrediction; i < traversalValues.size(); i++) {
        const auto &[iteration, time] = traversalValues[i];
        iterationMultTime += iteration * time;
        iterationSum += iteration;
        iterationSquareSum += iteration * iteration;
        timeSum += time;
      }

      // cast integer to decimal because this division contains small numbers which would cause precision lose
      const double iterationMeanValue =
          static_cast<double>(iterationSum) / static_cast<double>(_evidenceFirstPrediction);
      const auto timeMeanValue = timeSum / _evidenceFirstPrediction;

      // cast unsigned to signed because this difference can be negative
      const auto numerator = static_cast<long>(iterationMultTime) - static_cast<long>(iterationSum * timeMeanValue);
      const auto denominator =
          static_cast<double>(iterationSquareSum) - static_cast<double>(iterationSum * iterationMeanValue);

      // ((Sum iteration_i * time_i) - n * iterationMeanValue * timeMeanValue) / ((Sum iteration_i^2) - n *
      // iterationMeanValue ^ 2)
      const auto gradient = numerator / denominator;
      const auto yIntercept = timeMeanValue - gradient * iterationMeanValue;

      _configurationPredictions[configuration] = gradient * _iterationBeginTuningPhase + yIntercept;
      functionParams.clear();
      functionParams.emplace_back(gradient);
      functionParams.emplace_back(yIntercept);

      logPrediction(configuration);
    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = std::numeric_limits<long unsigned int>::max();
      _tooLongNotTestedSearchSpace.emplace(configuration);
      AutoPasLog(debug, "No traversal time prediction for {}", configuration.toString());
    }
  }
}

void PredictiveTuning::lagrangePolynomial() {
  for (const auto &configuration : _searchSpace) {
    const auto &traversalValues = _traversalTimesStorage[configuration];
    if (traversalValues.size() >= _evidenceFirstPrediction) {
      const auto lengthTTS = traversalValues.size() - 1;
      long prediction = 0;
      // calculates lagrange basis functions L_i and multiplies it with the tested time of p_i
      for (unsigned int i = 0; i < _evidenceFirstPrediction; i++) {
        long numerator = 1;
        long denominator = 1;
        const auto &[pointIIteration, pointITime] = traversalValues[lengthTTS - i];
        for (unsigned int j = 0; j < _evidenceFirstPrediction; j++) {
          if (i != j) {
            // cast unsigned to signed because this difference can be negative
            numerator *= static_cast<long>(_iterationBeginTuningPhase) - traversalValues[lengthTTS - j].first;
            denominator *= static_cast<long>(pointIIteration) - traversalValues[lengthTTS - j].first;
            // TODO: Look if this is sufficient or if there is a better solution
            // Should only be needed when the denominator overflows and therefore gets set to zero
            if (denominator == 0) {
              denominator = std::numeric_limits<long>::max();
              break;
            }
          }
        }
        prediction += numerator * pointITime / denominator;
      }
      _configurationPredictions[configuration] = prediction;
      logPrediction(configuration);
    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = std::numeric_limits<long unsigned int>::max();
      _tooLongNotTestedSearchSpace.emplace(configuration);
      AutoPasLog(debug, "No traversal time prediction for {}", configuration.toString());
    }
  }
}

void PredictiveTuning::newtonPolynomial() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    // if configuration was not tested in last tuning phase reuse prediction function.
    if ((_lastTest[configuration] != (_tuningPhaseCounter - 1)) and functionParams.size() == _evidenceFirstPrediction) {
      const auto &traversalValues = _traversalTimesStorage[configuration];
      auto lengthTTS = traversalValues.size() - 1;
      long prediction = 0;
      for (unsigned int i = 0; i < _evidenceFirstPrediction; i++) {
        auto interimValue = functionParams[i];
        for (unsigned int j = 0; j < i; j++) {
          // cast to double, because _configurationPredictionFunctionParams contains doubles
          interimValue *= static_cast<double>(_iterationBeginTuningPhase) -
                          traversalValues[lengthTTS - _evidenceFirstPrediction + j].first;
        }
        prediction += interimValue;
      }
      _configurationPredictions[configuration] = prediction;
      logPrediction(configuration);
    } else if (const auto &traversalValues = _traversalTimesStorage[configuration];
               traversalValues.size() >= _evidenceFirstPrediction) {
      std::vector<std::vector<double>> interimCalculation(_evidenceFirstPrediction);
      std::vector<size_t> iterationValues(_evidenceFirstPrediction);
      std::vector<double> coefficients(_evidenceFirstPrediction);
      const auto lengthTTS = traversalValues.size();
      auto lengthIthColumn = _evidenceFirstPrediction;

      // calculates a column of the newton interpolation
      for (unsigned int i = 0; i < _evidenceFirstPrediction; i++) {
        std::vector<double> ithColumn(lengthIthColumn);
        for (unsigned int j = 0; j < lengthIthColumn; j++) {
          if (i == 0) {
            auto &[iteration, runTime] = traversalValues[lengthTTS - _evidenceFirstPrediction + j];
            ithColumn[j] = runTime;
            iterationValues[j] = iteration;
          } else {
            // cast integer to decimal because this division contains small numbers which would cause precision lose
            ithColumn[j] = (interimCalculation[i - 1][j + 1] - interimCalculation[i - 1][j]) /
                           static_cast<double>(iterationValues[j + i] - iterationValues[j]);
          }
        }
        interimCalculation[i] = ithColumn;
        coefficients[i] = ithColumn.front();
        --lengthIthColumn;
      }

      size_t prediction = 0;
      for (unsigned int i = 0; i < _evidenceFirstPrediction; i++) {
        auto interimValue = coefficients[i];
        for (unsigned int j = 0; j < i; j++) {
          // cast to double, because coefficients contains doubles
          interimValue *= static_cast<double>(_iterationBeginTuningPhase - iterationValues[j]);
        }
        prediction += interimValue;
      }

      _configurationPredictions[configuration] = prediction;
      logPrediction(configuration);
      functionParams = coefficients;
    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = std::numeric_limits<long unsigned int>::max();
      _tooLongNotTestedSearchSpace.emplace(configuration);
      AutoPasLog(debug, "No traversal time prediction for {}", configuration.toString());
    }
  }
}

void PredictiveTuning::logPrediction(const Configuration &configuration) {
  // print config, prediction
  if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
    std::ostringstream ss;
    // print config
    ss << configuration.toString() << " : ";
    // print prediction
    ss << _configurationPredictions[configuration];
    AutoPasLog(debug, "Traversal time prediction for {}", ss.str());
  }
}

void PredictiveTuning::reselectOptimalSearchSpace() {
  // This is placed here, because there are no unnecessary removals
  for (const auto &configuration : _optimalSearchSpace) {
    _configurationPredictions.erase(configuration);
    _validSearchSpace.erase(configuration);
  }

  _optimalSearchSpace.clear();

  if (_validSearchSpace.size() == 1) {
    _optimalSearchSpace = _validSearchSpace;
    _currentConfig = _optimalSearchSpace.begin();
    return;
  }

  // checks if there are any configurations left that can be tested
  if (_configurationPredictions.empty()) {
    autopas::utils::ExceptionHandler::exception("Predictive Tuning: No valid configuration could be found");
  }

  const auto optimum = getOptimum(_configurationPredictions);

  if (_validSearchSpace.count(optimum->first) == 0) {
    autopas::utils::ExceptionHandler::exception("Predictive Tuning: No valid optimal configuration could be found");
  }

  _optimalSearchSpace.emplace(optimum->first);

  // selects configurations that are near the optimal prediction
  for (const auto &configuration : _validSearchSpace) {
    // Adds configurations that are within the _relativeOptimumRange
    if ((float)_configurationPredictions[configuration] / optimum->second <= _relativeOptimumRange) {
      _optimalSearchSpace.emplace(configuration);
      _tooLongNotTestedSearchSpace.erase(configuration);
    }
  }

  // sanity check
  if (_optimalSearchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("PredicitveTuning: No possible configuration prediction found!");
  }

  _currentConfig = _optimalSearchSpace.begin();
}

bool PredictiveTuning::tune(bool currentInvalid) {
  if (not currentInvalid) {
    _validConfigurationFound = true;
  }

  // repeat as long as traversals are not applicable or we run out of configs
  ++_currentConfig;

  if (_currentConfig == _searchSpace.end() or _currentConfig == _optimalSearchSpace.end()) {
    if (_validConfigurationFound) {
      if (_tooLongNotTestedSearchSpace.empty()) {
        selectOptimalConfiguration();
        _tuningPhaseCounter++;
        return false;
      } else {
        _currentConfig = _tooLongNotTestedSearchSpace.begin();
        return true;
      }
    } else {
      reselectOptimalSearchSpace();
    }
  } else if (_currentConfig == _tooLongNotTestedSearchSpace.end()) {
    selectOptimalConfiguration();
    _tuningPhaseCounter++;
    return false;
  }

  return true;
}

void PredictiveTuning::selectOptimalConfiguration() {
  if (_optimalSearchSpace.size() == 1) {
    _currentConfig = _optimalSearchSpace.begin();
    AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
    return;
  }
  if (_searchSpace.size() == 1) {
    _currentConfig = _searchSpace.begin();
    AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
    return;
  }

  // select the tested traversal times for the current tuning phase
  std::unordered_map<Configuration, size_t, ConfigHash> traversalTimes;
  // In the first couple iterations tune iterates through _searchSpace until predictions are made
  if (_optimalSearchSpace.empty()) {
    for (const auto &configuration : _searchSpace) {
      if (_lastTest[configuration] == _tuningPhaseCounter) {
        traversalTimes[configuration] = _traversalTimesStorage[configuration].back().second;
      }
    }
  } else {
    for (const auto &configuration : _optimalSearchSpace) {
      if (_lastTest[configuration] == _tuningPhaseCounter) {
        traversalTimes[configuration] = _traversalTimesStorage[configuration].back().second;
      }
    }
    for (const auto &configuration : _tooLongNotTestedSearchSpace) {
      if (_lastTest[configuration] == _tuningPhaseCounter) {
        traversalTimes[configuration] = _traversalTimesStorage[configuration].back().second;
      }
    }
  }

  // Time measure strategy
  if (traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "PredictiveTuning: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  const auto optimum = getOptimum(traversalTimes);

  _currentConfig = _searchSpace.find(optimum->first);

  // sanity check
  if (_currentConfig == _searchSpace.end() or _currentConfig == _optimalSearchSpace.end() or
      _currentConfig == _tooLongNotTestedSearchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning: Optimal configuration not found in list of configurations!");
  }

  if (_relativeBlacklistRange != 0 and not _notTestedYet.empty()) {
    blacklistBadConfigurations();
  }

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

void PredictiveTuning::blacklistBadConfigurations() {
  auto optimumTime = static_cast<double>(_traversalTimesStorage[getCurrentConfiguration()].back().second);
  auto configurationIter = _notTestedYet.begin();
  while (configurationIter != _notTestedYet.end()) {
    if (_lastTest[*configurationIter] == _tuningPhaseCounter) {
      if (_traversalTimesStorage[*configurationIter].back().second / optimumTime > _relativeBlacklistRange) {
        _searchSpace.erase(*configurationIter);
        _traversalTimesStorage.erase(*configurationIter);
        _lastTest.erase(*configurationIter);
      }
      configurationIter = _notTestedYet.erase(configurationIter);
    } else {
      configurationIter++;
    }
  }
}

}  // namespace autopas
