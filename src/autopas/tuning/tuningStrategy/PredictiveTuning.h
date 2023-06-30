/**
 * @file PredictiveTuning.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <limits>
#include <set>
#include <utility>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/logging/PredictionLogger.h"
#include "tuning/searchSpace/EvidenceCollection.h"

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
class PredictiveTuning final : public TuningStrategyInterface {
 public:
  /**
   * Shorthand for the type used to store predictions in.
   * Vector< Configuration, PredictionValue >
   */
  using PredictionsType = std::vector<std::tuple<Configuration, long>>;

  /**
   * Constructor for the PredictiveTuning that generates the search space from the allowed options.
   * @param relativeOptimum
   * @param maxTuningIterationsWithoutTest
   * @param testsUntilFirstPrediction
   * @param extrapolationMethodOption
   * @param outputSuffix
   */
  PredictiveTuning(const std::set<Configuration> &searchSpace, double relativeOptimum,
                   unsigned int maxTuningIterationsWithoutTest, unsigned int testsUntilFirstPrediction,
                   ExtrapolationMethodOption extrapolationMethodOption, const std::string &outputSuffix = "")
      : _validSearchSpace(searchSpace),
        _relativeOptimumRange(relativeOptimum),
        _maxTuningPhasesWithoutTest(maxTuningIterationsWithoutTest),
        _extrapolationMethod(extrapolationMethodOption),
        _minNumberOfEvidence(
            extrapolationMethodOption == ExtrapolationMethodOption::linePrediction ? 2 : testsUntilFirstPrediction),
        _predictionLogger(outputSuffix) {}

  /**
   * Constructor for the PredictiveTuning that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit PredictiveTuning(const std::set<Configuration> &allowedConfigurations) {}

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue, const EvidenceCollection &evidence) override;

  void rejectConfigurationIndefinitely(const Configuration &configuration) override {
    _validSearchSpace.erase(configuration);
    _tooLongNotTestedSearchSpace.erase(configuration);
  };

 private:
  /**
   * For all given configuration use the chosen extrapolation method to calculate a prediction based on the
   * given evidence.
   *
   * @param iteration The iteration for which to calculate the predictions.
   * @param tuningPhase The current tuning phase number.
   * @param configurations Configurations for which to calculate predictions.
   * @param evidenceCollection Data on which all extrapolation is based on.
   * @return For each configuration a predictions.
   */
  PredictiveTuning::PredictionsType calculatePredictions(size_t iteration, size_t tuningPhase,
                                                         const std::vector<Configuration> &configurations,
                                                         const autopas::EvidenceCollection &evidenceCollection);
  /**
   * Predicts the traversal time by placing a line through the last two traversal points and calculating the prediction
   * for the current time.
   */
  long linePrediction(size_t iteration, size_t tuningPhase, const Configuration &configuration,
                      const std::vector<Evidence> &evidenceVec);
  /**
   * Predicts the traversal time by creating a function that places a line through the data points and calculating the
   * prediction for the current time.
   */
  long linearRegression(size_t iteration, size_t tuningPhase, const Configuration &configuration,
                        const std::vector<Evidence> &evidenceVec);
  /**
   * Creates a polynomial function using Newton's method of finite differences and with this function the prediction is
   * calculated.
   */
  long newtonPolynomial(size_t iteration, size_t tuningPhase, const Configuration &configuration,
                        const std::vector<Evidence> &evidenceVec);

  /**
   * Error value used as a placeholder for the predictions of configurations that are not predicted.
   */
  constexpr static long _predictionErrorValue{std::numeric_limits<long>::max()};

  /**
   * Placeholder value used when a prediction overflows.
   */
  constexpr static long _predictionOverflowValue{std::numeric_limits<long>::max() - 1};

  /**
   * Placeholder value used when a prediction overflows.
   */
  constexpr static long _predictionUnderflowValue{1l};

  /**
   * A Map that for each configuration stores the function for the prediction to reuse it if no new traversal time was
   * added in the last tuning  phase. The way that function is stored depends on the prediction method, hence it is a
   * vector:
   * Line Prediction: Gradient and last evidence
   * Linear Regression: Gradient and iteration
   * Newton: Vector of coefficients
   */
  std::unordered_map<Configuration, std::vector<double>, ConfigHash> _predictionFunctionParameters{};
  /**
   * Contains the configuration that have not been tested for a period of time and are going to be tested.
   */
  std::set<Configuration> _tooLongNotTestedSearchSpace{};
  /**
   * Contains the configurations that are not marked invalid in the current tuning phase.
   */
  std::set<Configuration> _validSearchSpace;
  /**
   * Factor of the range of the optimal configurations for the optimalSearchSpace.
   */
  double _relativeOptimumRange{1.2};
  /**
   * If a config is not tested for this number of tuning phases test it again to make predictions more reliable.
   */
  unsigned int _maxTuningPhasesWithoutTest{100};

  /**
   * Stores the extrapolation method that is going to be used for the traversal time predictions.
   */
  ExtrapolationMethodOption _extrapolationMethod{ExtrapolationMethodOption::linearRegression};
  /**
   * Stores the number of tests that have to be made until the first prediction.
   * This number also determines how much evidence is used for the prediction and for a polynomial extrapolation method,
   * this is also the degree of the polynomial.
   */
  unsigned int _minNumberOfEvidence{3};

  PredictionLogger _predictionLogger;
};

}  // namespace autopas
