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
#include "autopas/utils/Math.h"
#include "autopas/utils/logging/PredictionLogger.h"

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
   * @param outputSuffix
   */
  PredictiveTuning(const std::set<ContainerOption> &allowedContainerOptions,
                   const std::set<double> &allowedCellSizeFactors,
                   const std::set<TraversalOption> &allowedTraversalOptions,
                   const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                   const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                   const std::set<Newton3Option> &allowedNewton3Options, double relativeOptimum,
                   unsigned int maxTuningIterationsWithoutTest, double relativeRangeForBlacklist,
                   unsigned int testsUntilFirstPrediction, ExtrapolationMethodOption extrapolationMethodOption,
                   const std::string &outputSuffix = "")
      : SetSearchSpaceBasedTuningStrategy(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                                          allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options),
        _relativeOptimumRange(relativeOptimum),
        _maxTuningPhasesWithoutTest(maxTuningIterationsWithoutTest),
        _relativeBlacklistRange(relativeRangeForBlacklist),
        _extrapolationMethod(extrapolationMethodOption),
        _evidenceFirstPrediction(
            extrapolationMethodOption == ExtrapolationMethodOption::linePrediction ? 2 : testsUntilFirstPrediction),
        _predictionLogger(outputSuffix) {
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
      : SetSearchSpaceBasedTuningStrategy(std::move(allowedConfigurations)) {}

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

  /**
   * Getter for predicted runtimes.
   * @return _configurationPredictions
   *
   * @note only used for unit tests.
   */
  const std::unordered_map<Configuration, long, ConfigHash> &getConfigurationPredictions() const {
    return _configurationPredictions;
  };

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
   * Creates a polynomial function using Newton's method of finite differences and with this function the prediction is
   * calculated.
   */
  inline void newtonPolynomial();
  /**
   * Selects a new search space based on previous observations.
   * This method needs only to be called if everything in the previous _optimalSearchSpace was invalid.
   * Invalid configurations are discarded and
   * Creates a new optimalSearchSpace if every configuration in the previous one was invalid.
   */
  inline void reselectOptimalSearchSpace();
  /**
   * Purge configurations from the search space based on their performance in this tuning phase compared to the optimum.
   */
  inline void blacklistBadConfigurations();

  /**
   * Error value used as a placeholder for the predictions of configurations that are not predicted.
   */
  constexpr static long _predictionErrorValue{std::numeric_limits<long>::max()};

  /**
   * Placeholder value used when a prediciton overflows.
   */
  constexpr static long _predictionOverflowValue{std::numeric_limits<long>::max() - 1};

  /**
   * Placeholder value used when a prediciton overflows.
   */
  constexpr static long _predictionUnderflowValue{1l};

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
  std::unordered_map<Configuration, long, ConfigHash> _configurationPredictions{};
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
   * Stores the last tuning phase in which configuration was tested.
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

  PredictionLogger _predictionLogger;
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
    case ExtrapolationMethodOption::newton: {
      newtonPolynomial();
      break;
    }
  }
  // if AutoPas is compiled without -DAUTOPAS_LOG_PREDICTIONS this does nothing
  _predictionLogger.logAllPredictions(_searchSpace, _configurationPredictions, _predictionErrorValue,
                                      _tuningPhaseCounter);
}

void PredictiveTuning::linePrediction() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    // if configuration was not tested in last tuning phase reuse prediction function.
    // also make sure prediction function is of the right type (=correct amount of params)
    if (_lastTest[configuration] != (_tuningPhaseCounter - 1) and functionParams.size() == 2) {
      const auto delta = _firstIterationOfTuningPhase - _traversalTimesStorage[configuration].back().first;

      // gradient * delta + last point
      _configurationPredictions[configuration] =
          static_cast<long>(functionParams[0] * static_cast<double>(delta) + functionParams[1]);
    } else
        // if there is enough evidence calculate new prediction function
        if (const auto &traversalValues = _traversalTimesStorage[configuration];
            traversalValues.size() >= _evidenceFirstPrediction) {
      const auto &[traversal1Iteration, traversal1Time] = traversalValues[traversalValues.size() - 1];
      const auto &[traversal2Iteration, traversal2Time] = traversalValues[traversalValues.size() - 2];

      const long gradient = static_cast<long>(traversal1Time - traversal2Time) /
                            static_cast<long>(traversal1Iteration - traversal2Iteration);
      const long delta = static_cast<long>(_firstIterationOfTuningPhase) - static_cast<long>(traversal1Iteration);

      const long change = utils::Math::safeMul(gradient, delta);

      // this might overflow so use safeAdd.
      const long newValue =
          utils::Math::safeAdd(traversal1Time, change, _predictionUnderflowValue, _predictionOverflowValue);
      // Do not accept values smaller zero.
      _configurationPredictions[configuration] = newValue < 0 ? _predictionUnderflowValue : newValue;

      functionParams.clear();
      functionParams.emplace_back(gradient);
      functionParams.emplace_back(traversal1Time);
    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = _predictionErrorValue;
      _tooLongNotTestedSearchSpace.emplace(configuration);
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
      _configurationPredictions[configuration] =
          static_cast<long>(functionParams[0] * static_cast<double>(_firstIterationOfTuningPhase) + functionParams[1]);
    } else if (const auto &traversalValues = _traversalTimesStorage[configuration];
               traversalValues.size() >= _evidenceFirstPrediction) {
      // we need signed types because calculation of the gradient might have negative result
      long iterationMultTime = 0, timeSum = 0;
      size_t iterationSum = 0, iterationSquareSum = 0;

      bool numericOverflow = false;
      for (auto i = traversalValues.size() - _evidenceFirstPrediction; i < traversalValues.size(); i++) {
        const auto &[iteration, time] = traversalValues[i];
        long iterationMultTimeI = utils::Math::safeMul(static_cast<long>(iteration), time);
        iterationMultTime = utils::Math::safeAdd(iterationMultTime, iterationMultTimeI);
        // if any of the safe operations overflow we can directly move to the next config
        if (iterationMultTime == std::numeric_limits<decltype(iterationMultTime)>::max()) {
          numericOverflow = true;
          break;
        }
        iterationSum += iteration;
        // this will overflow at iteration 3 037 000 499
        iterationSquareSum += iteration * iteration;
        timeSum += time;
      }

      // if there is an overflow the actual prediction will also overflow so abort and continue.
      if (numericOverflow) {
        _configurationPredictions[configuration] = _predictionOverflowValue;
        continue;
      }

      // cast integer to decimal because this division contains small numbers which would cause precision lose
      const double iterationMeanValue =
          static_cast<double>(iterationSum) / static_cast<double>(_evidenceFirstPrediction);
      const long timeMeanValue = timeSum / _evidenceFirstPrediction;

      const long numerator = iterationMultTime - static_cast<long>(iterationSum) * timeMeanValue;
      const double denominator =
          static_cast<double>(iterationSquareSum) - static_cast<double>(iterationSum) * iterationMeanValue;

      // ((Sum iteration_i * time_i) - n * iterationMeanValue * timeMeanValue) / ((Sum iteration_i^2) - n *
      // iterationMeanValue ^ 2)
      const auto gradient = static_cast<double>(numerator) / denominator;

      const auto change = static_cast<long>(gradient * (_firstIterationOfTuningPhase - iterationMeanValue));
      // check if prediction runs into over or underflow.
      const long newValue =
          utils::Math::safeAdd(change, timeMeanValue, _predictionUnderflowValue, _predictionOverflowValue);
      // Do not accept values smaller zero.
      _configurationPredictions[configuration] = newValue < 0 ? _predictionUnderflowValue : newValue;

      functionParams.clear();
      functionParams.emplace_back(gradient);
      const auto yIntercept = static_cast<double>(timeMeanValue) - gradient * iterationMeanValue;
      functionParams.emplace_back(yIntercept);

    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = _predictionErrorValue;
      _tooLongNotTestedSearchSpace.emplace(configuration);
    }
  }
}

void PredictiveTuning::newtonPolynomial() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    bool numericOverflow = false;

    // if configuration was not tested in last tuning phase reuse prediction function.
    if ((_lastTest[configuration] != (_tuningPhaseCounter - 1)) and functionParams.size() == _evidenceFirstPrediction) {
      const auto &traversalValues = _traversalTimesStorage[configuration];
      auto lengthTTS = traversalValues.size() - 1;
      long prediction = 0;

      for (unsigned int i = 0; i < _evidenceFirstPrediction; i++) {
        auto interimValue = functionParams[i];
        for (unsigned int j = 0; j < i; j++) {
          // cast to double, because _configurationPredictionFunctionParams contains doubles
          auto factor = static_cast<double>(_firstIterationOfTuningPhase) -
                        traversalValues[lengthTTS - _evidenceFirstPrediction + j].first;
          // check if interimValue *= factor will overflow
          if (interimValue > 0 and factor > std::numeric_limits<double>::max() / interimValue) {
            numericOverflow = true;
            break;
          }
          interimValue *= factor;
        }
        // check if prediction += interimValue will overflow
        if (numericOverflow or static_cast<decltype(prediction)>(interimValue) >
                                   (std::numeric_limits<decltype(prediction)>::max() - prediction)) {
          numericOverflow = true;
        }
        // if the inner loop overflowed no point in continuing the outer one.
        if (numericOverflow) {
          break;
        }

        prediction += static_cast<decltype(prediction)>(interimValue);
      }

      if (numericOverflow) {
        _configurationPredictions[configuration] = _predictionOverflowValue;
      } else if (prediction < 0) {
        _configurationPredictions[configuration] = _predictionUnderflowValue;
      } else {
        _configurationPredictions[configuration] = prediction;
      }
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
            const auto &[iteration, runTime] = traversalValues[lengthTTS - _evidenceFirstPrediction + j];
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

      long prediction = 0;
      for (unsigned int i = 0; i < _evidenceFirstPrediction; i++) {
        auto interimValue = coefficients[i];
        for (unsigned int j = 0; j < i; j++) {
          // cast to double, because coefficients contains doubles
          auto factor = static_cast<double>(_firstIterationOfTuningPhase - iterationValues[j]);
          // check if interimValue *= factor will overflow
          if (interimValue > 0 and factor > std::numeric_limits<double>::max() / interimValue) {
            numericOverflow = true;
            break;
          }
          interimValue *= factor;
        }
        // check if prediction += interimValue will overflow
        if (numericOverflow or interimValue > static_cast<double>(std::numeric_limits<decltype(prediction)>::max()) or
            static_cast<decltype(prediction)>(interimValue) >
                (std::numeric_limits<decltype(prediction)>::max() - prediction)) {
          numericOverflow = true;
        }
        // if the inner loop overflowed no point in continuing the outer one.
        if (numericOverflow) {
          break;
        }

        prediction += interimValue;
      }

      if (numericOverflow) {
        _configurationPredictions[configuration] = _predictionOverflowValue;
      } else if (prediction < 0) {
        _configurationPredictions[configuration] = _predictionUnderflowValue;
      } else {
        _configurationPredictions[configuration] = prediction;
      }

      functionParams = coefficients;
    } else {
      // When a configuration was not yet tested twice.
      _configurationPredictions[configuration] = _predictionErrorValue;
      _tooLongNotTestedSearchSpace.emplace(configuration);
    }
  }
}

void PredictiveTuning::reselectOptimalSearchSpace() {
  // since this function is called only when the whole _optimalSearchSpace is invalid
  // mark all of it as invalid and remove it from the predictions
  for (const auto &configuration : _optimalSearchSpace) {
    _configurationPredictions.erase(configuration);
    _validSearchSpace.erase(configuration);
  }

  _optimalSearchSpace.clear();

  // shortcut if only one valid configuration is left
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

  // get optimum of remaining predictions
  const auto &[optimalConfig, optimalRuntime] = *getOptimum(_configurationPredictions);

  // abort if optimum is invalid
  if (_validSearchSpace.count(optimalConfig) == 0) {
    autopas::utils::ExceptionHandler::exception(
        "PredictiveTuning::reselectOptimalSearchSpace() : No valid optimal configuration could be found");
  }

  // rebuild _optimalSearchSpace from valid configurations that have similar performance as the optimum
  _optimalSearchSpace.emplace(optimalConfig);
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
  } else if (_searchSpace.size() == 1) {
    _currentConfig = _searchSpace.begin();
  } else {
    // select the tested traversal times for the current tuning phase
    std::unordered_map<Configuration, long, ConfigHash> traversalTimes{};
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
    if (_relativeBlacklistRange != 0) {
      blacklistBadConfigurations();
    }
  }
}

void PredictiveTuning::blacklistBadConfigurations() {
  const auto optimumTime = _traversalTimesStorage[getCurrentConfiguration()].back().second;
  const auto blacklistThreshold = _relativeBlacklistRange * optimumTime;
  for (auto configurationIter = _searchSpace.begin(); configurationIter != _searchSpace.end();) {
    // blacklist if the configuration was tested in this phase and is above the threshold
    if (_lastTest[*configurationIter] == _tuningPhaseCounter and
        _traversalTimesStorage[*configurationIter].back().second > blacklistThreshold) {
      _traversalTimesStorage.erase(*configurationIter);
      _lastTest.erase(*configurationIter);
      configurationIter = _searchSpace.erase(configurationIter);
    } else {
      configurationIter++;
    }
  }
}

}  // namespace autopas
