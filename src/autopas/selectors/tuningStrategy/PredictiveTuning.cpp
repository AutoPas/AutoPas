/**
 * @file PredictiveTuning.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "PredictiveTuning.h"

void autopas::PredictiveTuning::selectOptimalSearchSpace() {
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

void autopas::PredictiveTuning::calculatePredictions() {
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

void autopas::PredictiveTuning::linePrediction() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    // if configuration was not tested in last tuning phase reuse prediction function.
    // also make sure prediction function is of the right type (=correct amount of params)
    if (_lastTest[configuration] != (_tuningPhaseCounter - 1) and functionParams.size() == 2) {
      const auto delta = _firstIterationOfTuningPhase - _traversalTimesStorage[configuration].back().first;

      // gradient * delta + last point
      const double newValue = utils::Math::safeMul(functionParams[0], static_cast<double>(delta)) + functionParams[1];
      // check if multiplication overflowed and then explicitly set error value. No cast to avoid rounding errors.
      _configurationPredictions[configuration] =
          newValue == std::numeric_limits<double>::max() ? _predictionOverflowValue : static_cast<long>(newValue);
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

void autopas::PredictiveTuning::linearRegression() {
  for (const auto &configuration : _searchSpace) {
    auto &functionParams = _configurationPredictionFunctionParams[configuration];
    // if configuration was not tested in last tuning phase reuse prediction function.
    // checks is _configurationPredictionFunctionParams has a linear function saved and only then it can reuse it
    if ((_lastTest[configuration] != (_tuningPhaseCounter - 1)) and functionParams.size() == 2) {
      // if configuration was not tested in last tuning phase reuse prediction function.
      // gradient * iteration + y-intercept
      const double newValue =
          utils::Math::safeMul(functionParams[0], static_cast<double>(_firstIterationOfTuningPhase)) +
          functionParams[1];
      // check if multiplication overflowed and then explicitly set error value. No cast to avoid rounding errors.
      _configurationPredictions[configuration] =
          newValue == std::numeric_limits<double>::max() ? _predictionOverflowValue : static_cast<long>(newValue);
    } else if (const auto &traversalValues = _traversalTimesStorage[configuration];
               traversalValues.size() >= _evidenceFirstPrediction) {
      // we need signed types because calculation of the gradient might have negative result
      long iterationMultTime = 0;
      long timeSum = 0;
      size_t iterationSum = 0;
      size_t iterationSquareSum = 0;

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

      const auto numerator = static_cast<double>(iterationMultTime - static_cast<long>(iterationSum) * timeMeanValue);
      const double denominator =
          static_cast<double>(iterationSquareSum) - static_cast<double>(iterationSum) * iterationMeanValue;

      // ((Sum iteration_i * time_i) - n * iterationMeanValue * timeMeanValue) / ((Sum iteration_i^2) - n *
      // iterationMeanValue ^ 2)
      const auto gradient = numerator / denominator;

      const auto change =
          static_cast<long>(utils::Math::safeMul(gradient, (_firstIterationOfTuningPhase - iterationMeanValue)));
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

void autopas::PredictiveTuning::newtonPolynomial() {
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

void autopas::PredictiveTuning::reselectOptimalSearchSpace() {
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

bool autopas::PredictiveTuning::tune(bool currentInvalid) {
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

void autopas::PredictiveTuning::selectOptimalConfiguration() {
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

void autopas::PredictiveTuning::blacklistBadConfigurations() {
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
