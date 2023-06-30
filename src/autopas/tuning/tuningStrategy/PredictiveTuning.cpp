/**
 * @file PredictiveTuning.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "PredictiveTuning.h"

#include <algorithm>
#include <vector>

#include "tuning/searchSpace/Evidence.h"
#include "tuning/searchSpace/EvidenceCollection.h"

namespace autopas {

PredictiveTuning::PredictiveTuning(double relativeOptimum,
                                   unsigned int maxTuningIterationsWithoutTest, unsigned int testsUntilFirstPrediction,
                                   ExtrapolationMethodOption extrapolationMethodOption, const std::string &outputSuffix)
    :
      _relativeOptimumRange(relativeOptimum),
      _maxTuningPhasesWithoutTest(maxTuningIterationsWithoutTest),
      _extrapolationMethod(extrapolationMethodOption),
      _minNumberOfEvidence(
          extrapolationMethodOption == ExtrapolationMethodOption::linePrediction ? 2 : testsUntilFirstPrediction),
      _predictionLogger(outputSuffix) {}

PredictiveTuning::PredictionsType PredictiveTuning::calculatePredictions(
    size_t iteration, size_t tuningPhase, const std::vector<Configuration> &configurations,
    const EvidenceCollection &evidenceCollection) {
  PredictiveTuning::PredictionsType predictions;
  predictions.reserve(configurations.size());

  for (const auto &configuration : configurations) {
    const auto &evidenceVec = evidenceCollection.getEvidence(configuration);
    const auto predictionValue = [&]() {
      switch (_extrapolationMethod) {
        case ExtrapolationMethodOption::linePrediction: {
          return linePrediction(iteration, tuningPhase, configuration, evidenceVec);
        }
        case ExtrapolationMethodOption::linearRegression: {
          return linearRegression(iteration, tuningPhase, configuration, evidenceVec);
        }
        case ExtrapolationMethodOption::newton: {
          return newtonPolynomial(iteration, tuningPhase, configuration, evidenceVec);
        }
      }
    }();
    predictions.emplace_back(configuration, predictionValue);
  }
  // if AutoPas is compiled without -DAUTOPAS_LOG_PREDICTIONS this does nothing
  _predictionLogger.logAllPredictions(predictions, _predictionErrorValue, tuningPhase);

  return predictions;
}

long PredictiveTuning::linePrediction(size_t iteration, size_t tuningPhase, const Configuration &configuration,
                                      const std::vector<Evidence> &evidenceVec) {
  auto &functionParams = _predictionFunctionParameters[configuration];
  // if configuration was not tested in last tuning phase reuse prediction function.
  // also make sure prediction function is of the right type (=correct amount of params)
  if (evidenceVec.back().tuningPhase != (tuningPhase - 1) and functionParams.size() == 2) {
    const auto delta = iteration - evidenceVec.back().iteration;
    // gradient * delta + last point
    const double prediction = utils::Math::safeMul(functionParams[0], static_cast<double>(delta)) + functionParams[1];
    // check if multiplication overflowed and then explicitly set error value. No cast to avoid rounding errors.
    const auto predictionOverflowChecked =
        prediction == std::numeric_limits<double>::max() ? _predictionOverflowValue : static_cast<long>(prediction);
    return predictionOverflowChecked;

  } else
      // if there is enough evidenceVec calculate new prediction function
      if (evidenceVec.size() >= _minNumberOfEvidence) {
    const auto &[traversal1Iteration, traversal1TuningPhase, traversal1Time] = evidenceVec[evidenceVec.size() - 1];
    const auto &[traversal2Iteration, traversal2TuningPhase, traversal2Time] = evidenceVec[evidenceVec.size() - 2];

    const long gradient = static_cast<long>(traversal1Time - traversal2Time) /
                          static_cast<long>(traversal1Iteration - traversal2Iteration);
    const long delta = static_cast<long>(iteration) - static_cast<long>(traversal1Iteration);

    const long change = utils::Math::safeMul(gradient, delta);

    // this might overflow so use safeAdd.
    const long prediction =
        utils::Math::safeAdd(traversal1Time, change, _predictionUnderflowValue, _predictionOverflowValue);
    // Do not accept values smaller zero.
    const auto predictionUnderflowChecked = prediction < 0 ? _predictionUnderflowValue : prediction;

    // update _predictionFunctionParameters
    functionParams.clear();
    functionParams.emplace_back(gradient);
    functionParams.emplace_back(traversal1Time);

    return predictionUnderflowChecked;
  } else {
    // When a configuration was not yet measured enough to make a prediction insert a placeholder.
    _tooLongNotTestedSearchSpace.emplace(configuration);
    return _predictionErrorValue;
  }
}

long PredictiveTuning::linearRegression(size_t iteration, size_t tuningPhase, const Configuration &configuration,
                                        const std::vector<Evidence> &evidenceVec) {
  auto &functionParams = _predictionFunctionParameters[configuration];
  // if configuration was not tested in last tuning phase reuse prediction function.
  // checks is _configurationPredictionFunctionParams has a linear function saved and only then it can reuse it
  if ((evidenceVec.back().tuningPhase != (tuningPhase - 1)) and functionParams.size() == 2) {
    // if configuration was not tested in last tuning phase reuse prediction function.
    // gradient * iteration + y-intercept
    const double prediction =
        utils::Math::safeMul(functionParams[0], static_cast<double>(iteration)) + functionParams[1];
    // check if multiplication overflowed and then explicitly set error value. No cast to avoid rounding errors.
    const auto predictionOverflowChecked =
        prediction == std::numeric_limits<double>::max() ? _predictionOverflowValue : static_cast<long>(prediction);
    return predictionOverflowChecked;

  } else
      // if there is enough evidenceVec calculate new prediction function
      if (evidenceVec.size() >= _minNumberOfEvidence) {
    // we need signed types because calculation of the gradient might have negative result
    long iterationMultTime = 0;
    long timeSum = 0;
    size_t iterationSum = 0;
    size_t iterationSquareSum = 0;

    bool numericOverflow = false;
    for (auto i = evidenceVec.size() - _minNumberOfEvidence; i < evidenceVec.size(); i++) {
      const auto &[evidenceIteration, evidenceTuningPhase, evidenceValue] = evidenceVec[i];
      const auto iterationMultTimeI = utils::Math::safeMul(static_cast<long>(evidenceIteration), evidenceValue);
      iterationMultTime = utils::Math::safeAdd(iterationMultTime, iterationMultTimeI);
      // if any of the safe operations overflow we can directly move to the next config
      if (iterationMultTime == std::numeric_limits<decltype(iterationMultTime)>::max()) {
        numericOverflow = true;
        break;
      }
      iterationSum += evidenceIteration;
      iterationSquareSum += evidenceIteration * evidenceIteration;
      timeSum += evidenceValue;
    }

    // if there is an overflow the actual prediction will also overflow so abort and continue.
    if (numericOverflow) {
      return _predictionOverflowValue;
    }

    // cast integer to decimal because this division contains small numbers which would cause precision lose
    const double iterationMeanValue = static_cast<double>(iterationSum) / static_cast<double>(_minNumberOfEvidence);
    const long timeMeanValue = timeSum / _minNumberOfEvidence;

    const auto numerator = static_cast<double>(iterationMultTime - static_cast<long>(iterationSum) * timeMeanValue);
    const double denominator =
        static_cast<double>(iterationSquareSum) - static_cast<double>(iterationSum) * iterationMeanValue;

    // ((Sum iteration_i * time_i) - n * iterationMeanValue * timeMeanValue) / ((Sum iteration_i^2) - n *
    // iterationMeanValue ^ 2)
    const auto gradient = numerator / denominator;

    const auto change =
        static_cast<long>(utils::Math::safeMul(gradient, (static_cast<double>(iteration) - iterationMeanValue)));
    // check if prediction runs into over or underflow.
    const long prediction =
        utils::Math::safeAdd(change, timeMeanValue, _predictionUnderflowValue, _predictionOverflowValue);
    // Do not accept values smaller zero.
    const auto predictionUnderflowChecked = prediction < 0 ? _predictionUnderflowValue : prediction;

    functionParams.clear();
    functionParams.emplace_back(gradient);
    const auto yIntercept = static_cast<double>(timeMeanValue) - gradient * iterationMeanValue;
    functionParams.emplace_back(yIntercept);

    return predictionUnderflowChecked;
  } else {
    // When a configuration was not yet measured enough to make a prediction insert a placeholder.
    _tooLongNotTestedSearchSpace.emplace(configuration);
    return _predictionErrorValue;
  }
}

long PredictiveTuning::newtonPolynomial(size_t iteration, size_t tuningPhase, const Configuration &configuration,
                                        const std::vector<Evidence> &evidenceVec) {
  auto &functionParams = _predictionFunctionParameters[configuration];
  bool numericOverflow = false;

  // if configuration was not tested in last tuning phase reuse prediction function.
  if ((evidenceVec.back().tuningPhase != (tuningPhase - 1)) and functionParams.size() == _minNumberOfEvidence) {
    const auto numberOfEvidenceMinusOne = evidenceVec.size() - 1;
    long prediction = 0;

    for (unsigned int i = 0; i < _minNumberOfEvidence; i++) {
      auto interimValue = functionParams[i];
      for (unsigned int j = 0; j < i; j++) {
        // cast to double, because _configurationPredictionFunctionParams contains doubles
        const auto factor =
            static_cast<double>(iteration - evidenceVec[numberOfEvidenceMinusOne - _minNumberOfEvidence + j].iteration);
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
      return _predictionOverflowValue;
    } else if (prediction < 0) {
      return _predictionUnderflowValue;
    } else {
      return prediction;
    }

  } else if (evidenceVec.size() >= _minNumberOfEvidence) {
    std::vector<std::vector<double>> interimCalculation(_minNumberOfEvidence);
    std::vector<size_t> iterationValues(_minNumberOfEvidence);
    std::vector<double> coefficients(_minNumberOfEvidence);
    const auto numberOfEvidence = evidenceVec.size();
    auto lengthIthColumn = _minNumberOfEvidence;

    // calculates a column of the newton interpolation
    for (unsigned int i = 0; i < _minNumberOfEvidence; i++) {
      std::vector<double> ithColumn(lengthIthColumn);
      for (unsigned int j = 0; j < lengthIthColumn; j++) {
        if (i == 0) {
          const auto &[evidenceIteration, evidenceTuningPhase, evidenceValue] =
              evidenceVec[numberOfEvidence - _minNumberOfEvidence + j];
          ithColumn[j] = static_cast<double>(evidenceValue);
          iterationValues[j] = evidenceIteration;
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
    for (unsigned int i = 0; i < _minNumberOfEvidence; i++) {
      auto interimValue = coefficients[i];
      for (unsigned int j = 0; j < i; j++) {
        // cast to double, because coefficients contains doubles
        auto factor = static_cast<double>(iteration - iterationValues[j]);
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

    functionParams = coefficients;

    if (numericOverflow) {
      return _predictionOverflowValue;
    } else if (prediction < 0) {
      return _predictionUnderflowValue;
    } else {
      return prediction;
    }

  } else {
    // When a configuration was not yet tested twice.
    _tooLongNotTestedSearchSpace.emplace(configuration);
    return _predictionErrorValue;
  }
}

void PredictiveTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                           const EvidenceCollection &evidence) {}

void PredictiveTuning::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                             const EvidenceCollection &evidenceCollection) {
  // collect all configurations that were not tested for too long
  for (const auto &conf : configQueue) {
    const auto numPhasesWithoutTest = tuningPhase - evidenceCollection.getEvidence(conf).back().tuningPhase;
    if (numPhasesWithoutTest >= _maxTuningPhasesWithoutTest) {
      _tooLongNotTestedSearchSpace.emplace(conf);
    }
  }

  const auto predictions = calculatePredictions(iteration, tuningPhase, configQueue, evidenceCollection);
  // find the best prediction
  const auto &[bestConf, bestPrediction] =
      *std::min_element(predictions.begin(), predictions.end(), [&](const auto &tupleA, const auto &tupleB) {
        const auto &[confA, predictionA] = tupleA;
        const auto &[confB, predictionB] = tupleB;
        return predictionA < predictionB;
      });
  const auto bestAcceptablePredictionValue =
      static_cast<long>(_relativeOptimumRange * static_cast<double>(bestPrediction));
  // Flush configurations and only re-insert those with good predictions.
  configQueue.clear();
  for (const auto &[conf, prediction] : predictions) {
    if (prediction < bestAcceptablePredictionValue) {
      configQueue.push_back(conf);
    }
  }
  // Also insert configurations that were not tested for too long and we want to check again
  for (const auto &conf : _tooLongNotTestedSearchSpace) {
    // don't create duplicates
    if (std::find(configQueue.begin(), configQueue.end(), conf) != configQueue.end()) {
      configQueue.push_back(conf);
    }
  }

  // sort configurations to minimize container conversion overhead.
  std::sort(configQueue.begin(), configQueue.end());
}

void PredictiveTuning::addEvidence(const Configuration &configuration, const Evidence & /*evidence*/) {
  // if this configuration was on the watchlist of configurations that were not tested for a while take it off.
  if (const auto confIter = _tooLongNotTestedSearchSpace.find(configuration);
      confIter != _tooLongNotTestedSearchSpace.end()) {
    _tooLongNotTestedSearchSpace.erase(confIter);
  }
}

void PredictiveTuning::rejectConfigurationIndefinitely(const Configuration &configuration) {
  _tooLongNotTestedSearchSpace.erase(configuration);
}
}  // namespace autopas