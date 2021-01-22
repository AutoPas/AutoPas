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
 * In every tuning phase, this strategy makes a prediction about the run time for every configuration.
 * Then only those are tested which have the best predictions. In the end, the configuration
 * that performed best during testing is selected.
 *
 * Predictions about run time are extrapolations from previous measurements. There are multiple
 * extrapolation methods available. Depending on this choice a certain number of tuning phases is
 * necessary where everything is tested.
 *
 * Additional features:
 * - Retesting of old configurations: Configurations are guaranteed to be retested after a given amount
 *   of tuning phases so that their predictions do not rely on data points that are too old.
 * - Blacklisting: Configurations that perform extremely bad will be completely disregarded forever.
 *
 * The strategy works by having multiple sets of configurations (e.g. the whole search space, optimal
 * search space, search space of configurations that were not tested for a long time). _currentConfig is
 * a iterator to any of them and might be switched between the sets depending on what is currently tested.
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
        _maxTuningPhasesWithoutTest(maxTuningIterationsWithoutTest),
        _relativeBlacklistRange(relativeRangeForBlacklist),
        _extrapolationMethod(extrapolationMethodOption),
        _evidenceFirstPrediction(
            extrapolationMethodOption == ExtrapolationMethodOption::linePrediction ? 2 : testsUntilFirstPrediction) {
    // sets traversalTimesStorage
    _traversalTimesStorage.reserve(_searchSpace.size());
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

  [[nodiscard]] inline long getEvidence(Configuration configuration) const override {
    // compute the average of times for this configuration
    auto times = _traversalTimesStorage.at(configuration);
    long result = 0;
    for (auto time : times) {
      result += time.second;
    }
    return result / times.size();
  }

  [[nodiscard]] inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void reset(size_t iteration) override {
    _configurationPredictions.clear();
    _optimalSearchSpace.clear();
    _tooLongNotTestedSearchSpace.clear();
    _validSearchSpace = _searchSpace;
    _validConfigurationFound = false;
    _firstIterationOfTuningPhase = iteration;

    selectOptimalSearchSpace();
  }

  inline bool tune(bool currentInvalid = false) override;

 private:
  /**
   * Selects the optimal (=fastest) configuration
   */
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
   * Selects a new search space based on previous observations.
   * Invalid configurations are discarded and
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
  inline void logPrediction(const Configuration &configuration);

  /**
   * For each configuration stores how long it took in which iteration.
   * Configuration -> vector<(iteration, runTime)>
   */
  std::unordered_map<Configuration, std::vector<std::pair<size_t, long>>, ConfigHash> _traversalTimesStorage{};
  /**
   * A Map that for each configuration stores the function for the prediction to reuse it if no new traversal time was
   * added in the last tuning  phase. The way that function is stored depends on the prediction method, hence it is a
   * vector:
   * Line Prediction: Gradient and last evidence
   * Linear Regression: Gradient and iteration
   * Newton: Vector of coefficients
   */
  std::unordered_map<Configuration, std::vector<double>, ConfigHash> _configurationPredictionFunctionParams{};
  /**
   * Mapping of configuration to their predicted run times.
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _configurationPredictions{};
  /**
   * Contains the configuration that are predicted to be potentially optimal and going to be tested.
   */
  std::set<Configuration> _optimalSearchSpace{};
  /**
   * Contains the configuration that have not been tested for a period of time and are going to be tested.
   */
  std::set<Configuration> _tooLongNotTestedSearchSpace{};
  /**
   * Contains the configurations that are not marked invalid in the current tuning phase.
   */
  std::set<Configuration> _validSearchSpace{};
  /**
   * Contains the configurations that have not been tested yet.
   */
  std::set<Configuration> _notTestedYet{};
  /**
   * Stores the the last tuning phase in which configuration was tested.
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _lastTest{};
  /**
   * Gets incremented after every completed tuning phase.
   * This is not zero-based, because when comparing with a value in _lastTest we would not be able to
   * distinguish easily if a configuration was tested in phase zero or not at all.
   */
  unsigned int _tuningPhaseCounter{1};
  /**
   * Stores the iteration at the beginning of a tuning phase.
   * Mainly used for the traversal time storage.
   */
  unsigned int _firstIterationOfTuningPhase{0};
  /**
   * Indicates if a valid configuration was found in the _optimalSearchSpace.
   */
  bool _validConfigurationFound{false};
  /**
   * Factor of the range of the optimal configurations for the optimalSearchSpace.
   */
  const double _relativeOptimumRange{1.2};
  /**
   * After not being tested this number of tuningPhases a configuration is being emplaced in _optimalSearchSpace.
   */
  const unsigned int _maxTuningPhasesWithoutTest{5};
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

void PredictiveTuning::selectOptimalSearchSpace() {
  // Special case: Search space is trivial or more evidence is needed before predictions can be made.
  if (_searchSpace.size() == 1 or _tuningPhaseCounter < _evidenceFirstPrediction + 1) {
    if (_searchSpace.empty()) {
      autopas::utils::ExceptionHandler::exception(
          "PredictiveTuning::selectOptimalSearchSpace() : No possible configuration prediction found!");
    }
    _currentConfig = _searchSpace.begin();
    return;
  }

  calculatePredictions();

  const auto &[optimalConfig, optimalRuntime] = *getOptimum(_configurationPredictions);
  const auto absoluteOptimalRange = _relativeOptimumRange * optimalRuntime;

  _optimalSearchSpace.emplace(optimalConfig);

  // build _optimalSearchSpace from all configurations that have similar performance as the optimum
  for (const auto &configuration : _searchSpace) {
    // Adds configurations that are within the _relativeOptimumRange
    if (static_cast<double>(_configurationPredictions[configuration]) <= absoluteOptimalRange) {
      _optimalSearchSpace.emplace(configuration);
      // ... or have not been tested for a long time
    } else if (_tuningPhaseCounter - _lastTest[configuration] > _maxTuningPhasesWithoutTest) {
      _tooLongNotTestedSearchSpace.emplace(configuration);
    }
  }

  // sanity check
  if (_optimalSearchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning::selectOptimalSearchSpace() : No possible configuration prediction found!");
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
    // also make sure prediction function is of the right type (=correct amount of params)
    if (_lastTest[configuration] != (_tuningPhaseCounter - 1) and functionParams.size() == 2) {
      const auto delta = _firstIterationOfTuningPhase - _traversalTimesStorage[configuration].back().first;

      // gradient * delta + last point
      _configurationPredictions[configuration] = functionParams[0] * delta + functionParams[1];
      logPrediction(configuration);
    } else
        // if there is enough evidence calculate new prediction function
        if (const auto &traversalValues = _traversalTimesStorage[configuration];
            traversalValues.size() >= _evidenceFirstPrediction) {
      const auto &[traversal1Iteration, traversal1Time] = traversalValues[traversalValues.size() - 1];
      const auto &[traversal2Iteration, traversal2Time] = traversalValues[traversalValues.size() - 2];

      const auto gradient = (traversal1Time - traversal2Time) / (traversal1Iteration - traversal2Iteration);
      const auto delta = _firstIterationOfTuningPhase - traversal1Iteration;

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
      _configurationPredictions[configuration] = functionParams[0] * _firstIterationOfTuningPhase + functionParams[1];
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

      _configurationPredictions[configuration] = gradient * _firstIterationOfTuningPhase + yIntercept;
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
            // cast unsigned to signed before subtraction, because this difference can be negative
            numerator *= static_cast<long>(_firstIterationOfTuningPhase) - traversalValues[lengthTTS - j].first;
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
          interimValue *= static_cast<double>(_firstIterationOfTuningPhase) -
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
          interimValue *= static_cast<double>(_firstIterationOfTuningPhase - iterationValues[j]);
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
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning::reselectOptimalSearchSpace() : No valid configuration could be found");
  }

  // optimumIter: first = config ; second = runTime
  const auto &[optimalConfig, optimalRuntime] = *getOptimum(_configurationPredictions);

  // abort if optimum is invalid
  if (_validSearchSpace.count(optimalConfig) == 0) {
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning::reselectOptimalSearchSpace() : No valid optimal configuration could be found");
  }

  _optimalSearchSpace.emplace(optimalConfig);

  // rebuild _optimalSearchSpace from valid configurations that have similar performance as the optimum
  const auto absoluteOptimumRange = _relativeOptimumRange * optimalRuntime;
  for (const auto &configuration : _validSearchSpace) {
    // Adds configurations that are within the _relativeOptimumRange
    if (_configurationPredictions[configuration] <= absoluteOptimumRange) {
      _optimalSearchSpace.emplace(configuration);
      _tooLongNotTestedSearchSpace.erase(configuration);
    }
  }

  // sanity check
  if (_optimalSearchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning::reselectOptimalSearchSpace() : No possible configuration prediction found!");
  }

  _currentConfig = _optimalSearchSpace.begin();
}

bool PredictiveTuning::tune(bool currentInvalid) {
  if (not currentInvalid) {
    _validConfigurationFound = true;
  }

  ++_currentConfig;

  // if we tested every config that needed testing (and are not iterating through _tooLongNotTestedSearchSpace)...
  if (_currentConfig == _searchSpace.end() or _currentConfig == _optimalSearchSpace.end()) {
    // ... and at least one config was valid
    if (_validConfigurationFound) {
      // ... and there is nothing that needs re-testing
      if (_tooLongNotTestedSearchSpace.empty()) {
        // finalize this tuning phase and settle on a configuration
        selectOptimalConfiguration();
        ++_tuningPhaseCounter;
        return false;
      } else {
        _currentConfig = _tooLongNotTestedSearchSpace.begin();
        return true;
      }
    } else {
      reselectOptimalSearchSpace();
    }
  } else if (_currentConfig == _tooLongNotTestedSearchSpace.end()) {
    // if we reach this point, we really tested everything there is to test
    selectOptimalConfiguration();
    ++_tuningPhaseCounter;
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
  std::unordered_map<Configuration, size_t, ConfigHash> traversalTimes{};
  // as long as there is too few evidence there is no optimal search space
  if (_optimalSearchSpace.empty()) {
    traversalTimes.reserve(_searchSpace.size());
    for (const auto &configuration : _searchSpace) {
      // check if this config was tests in the current phase
      if (_lastTest[configuration] == _tuningPhaseCounter) {
        traversalTimes[configuration] = _traversalTimesStorage[configuration].back().second;
      }
    }
  } else {
    traversalTimes.reserve(_optimalSearchSpace.size() + _tooLongNotTestedSearchSpace.size());
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
        "PredictiveTuning::selectOptimalConfiguration() : Trying to determine fastest configuration without any "
        "measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  const auto &[optimalConfig, _] = *getOptimum(traversalTimes);

  _currentConfig = _searchSpace.find(optimalConfig);

  // sanity check
  if (_currentConfig == _searchSpace.end() or _currentConfig == _optimalSearchSpace.end() or
      _currentConfig == _tooLongNotTestedSearchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning::selectOptimalConfiguration() : Optimal configuration not found in any set of "
        "configurations!");
  }

  // if a blacklist range is set blacklist the really bad stuff
  if (_relativeBlacklistRange != 0 and not _notTestedYet.empty()) {
    blacklistBadConfigurations();
  }

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

void PredictiveTuning::blacklistBadConfigurations() {
  const auto optimumTime = _traversalTimesStorage[getCurrentConfiguration()].back().second;
  const auto blacklistThreshold = _relativeBlacklistRange * optimumTime;
  for (auto configurationIter = _notTestedYet.begin(); configurationIter != _notTestedYet.end();) {
    if (_lastTest[*configurationIter] == _tuningPhaseCounter) {
      if (_traversalTimesStorage[*configurationIter].back().second > blacklistThreshold) {
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
