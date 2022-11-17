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

}  // namespace autopas
