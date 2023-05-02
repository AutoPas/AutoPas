/**
 * @file BayesianClusterSearch.h
 * @author Jan Nguyen
 * @date 12.04.20
 */

#pragma once

#include <limits>
#include <map>
#include <set>
#include <unordered_set>

#include "FullSearch.h"
#include "GaussianModel/GaussianCluster.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/FeatureVectorEncoder.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time while fixing discrete variables corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given
 * feature next.
 */
class BayesianClusterSearch : public TuningStrategyInterface {
  /**
   * The maximum number of attempts to sample an optimum.
   */
  constexpr static size_t maxAttempts = 10;
  /**
   * The factor for conversion from seconds to microseconds.
   */
  constexpr static double secondsPerMicroseconds = 1. / 1000000.;
  /**
   * Dimensions of the continuous tuples.
   * FeatureVector continuous dimensions + 1 (time)
   */
  constexpr static size_t continuousDims = FeatureVectorEncoder::tunableContinuousDims + 1;
  /**
   * Fixed noise
   */
  constexpr static double sigma = 0.01;

  /**
   * Iteration numbers get scaled down by a multiple of max evidence to ensure not all kernels become 0.
   */
  constexpr static double iterationScalePerMaxEvidence = 5.;

  /**
   * Distance between evidence is expected not to exceed a upper threshold.
   * We chose the threshold 7 because: exp(-r) < 0.1% for r > 7 but this
   * can be increased if necessary.
   */
  constexpr static double suggestedMaxDistance = 7.;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedVerletRebuildFrequencies
   * @param maxEvidence Stop tuning after given number of evidence provided.
   * @param predAcqFunction Acquisition function used for prediction while tuning.
   * @param outputSuffix Suffix for output logger.
   * @param predNumLHSamples Number of latin-hypercube-samples used to find a evidence with high predicted acquisition
   * @param seed Seed of random number generator (should only be used for tests)
   */
  explicit BayesianClusterSearch(
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(),
      const std::set<int> &allowedVerletRebuildFrequencies = std::set<int>({12, 24, 48}), size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      const std::string &outputSuffix = "", size_t predNumLHSamples = 50, unsigned long seed = std::random_device()())
      : _containerOptionsSet(allowedContainerOptions),
        _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
        _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _verletRebuildFrequencies((NumberSetFinite<int>(allowedVerletRebuildFrequencies)).clone()),
        _encoder(),
        _currentConfig(),
        _invalidConfigs(),
        _traversalTimes(),
        _rng(seed),
        _gaussianCluster({}, continuousDims, GaussianCluster::WeightFunction::evidenceMatchingScaledProbabilityGM,
                         sigma, _rng, GaussianCluster::defaultVecToString, outputSuffix),
        _neighbourFun([this](const Eigen::VectorXi &target) -> std::vector<std::pair<Eigen::VectorXi, double>> {
          return _encoder.clusterNeighboursManhattan1Container(target);
        }),
        _maxEvidence(maxEvidence),
        _predAcqFunction(predAcqFunction),
        _predNumLHSamples(predNumLHSamples),
        _firstTuningPhase(true),
        _currentIteration(0),
        _iterationScale(1. / (maxEvidence * iterationScalePerMaxEvidence)),
        _currentNumEvidence(0),
        _currentOptimalTime(std::numeric_limits<long>::max()),
        _fullSearch(allowedContainerOptions, {allowedCellSizeFactors.getMedian()}, allowedTraversalOptions,
                    allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options,
                    allowedVerletRebuildFrequencies) {
    if (predNumLHSamples <= 0) {
      utils::ExceptionHandler::exception(
          "BayesianSearch: Number of samples used for predictions must be greater than 0!");
    }

    for (const auto &containerOption : allowedContainerOptions) {
      // get all traversals of the container and restrict them to the allowed ones
      const std::set<TraversalOption> &allContainerTraversals =
          compatibleTraversals::allCompatibleTraversals(containerOption);
      std::set<TraversalOption> allowedAndApplicable;
      std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                            allContainerTraversals.begin(), allContainerTraversals.end(),
                            std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

      for (const auto &traversalOption : allowedAndApplicable) {
        // if load estimators are not applicable LoadEstimatorOption::none is returned.
        const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
            loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
        for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
          _containerTraversalEstimatorOptions.emplace_back(containerOption, traversalOption, loadEstimatorOption);
        }
      }
    }

    if (searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("BayesianClusterSearch: No valid configurations could be created.");
    }

    updateOptions();
    _gaussianCluster.setVectorToStringFun(
        [this](const GaussianModelTypes::VectorPairDiscreteContinuous &vec) -> std::string {
          return _encoder.convertFromCluster(vec).toString();
        });

    _currentConfig = _fullSearch.getCurrentConfiguration();
  }

  inline const Configuration &getCurrentConfiguration() const override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t iteration) override {
    if (_firstTuningPhase) {
      _fullSearch.addEvidence(time, iteration);
    }

    // store optimal evidence
    if (time < _currentOptimalTime) {
      _currentOptimalTime = time;
      _currentOptimalConfig = _currentConfig;
    }

    // encoded vector
    auto vec = _encoder.convertToCluster(_currentConfig, iteration * _iterationScale);
    // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
    // represent a maximization problem
    _gaussianCluster.addEvidence(vec, -time * secondsPerMicroseconds);

    _currentIteration = iteration;
    ++_currentNumEvidence;
    _currentAcquisitions.clear();

    _traversalTimes[_currentConfig] = time;
  }

  inline long getEvidence(Configuration configuration) const override { return _traversalTimes.at(configuration); }

  inline void reset(size_t iteration) override {
    size_t iterationSinceLastEvidence = iteration - _currentIteration;
    if (iterationSinceLastEvidence * _iterationScale > suggestedMaxDistance) {
      AutoPasLog(warn, "BayesianClusterSearch: Time since the last evidence may be too long.");
    }

    _currentIteration = iteration;
    _currentNumEvidence = 0;
    _currentOptimalTime = std::numeric_limits<long>::max();
    _currentAcquisitions.clear();

    if (_firstTuningPhase) {
      _fullSearch.reset(iteration);
      _currentConfig = _fullSearch.getCurrentConfiguration();
    } else {
      tune();
    }
  }

  inline bool tune(bool currentInvalid = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptionsSet; }

  inline bool searchSpaceIsTrivial() const override;

  inline bool searchSpaceIsEmpty() const override;

  inline bool smoothedHomogeneityAndMaxDensityNeeded() const override { return false; }

 private:
  /**
   * Generate n samples and predict their corresponding
   * acquisition function. The result is sorted depending on
   * the given acquisition function.
   * @param n numSamples
   * @param af acquisition function
   */
  inline void sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  /**
   * If allowed options are changed this functions should be called
   * to update encoder and clusters.
   */
  inline void updateOptions();

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
  std::unique_ptr<NumberSet<int>> _verletRebuildFrequencies;
  FeatureVectorEncoder _encoder;

  FeatureVector _currentConfig;
  /**
   * Currently sampled vectors and corresponding acquisition values.
   */
  std::vector<GaussianModelTypes::VectorPairDiscreteContinuous> _currentAcquisitions;
  /**
   * Configurations marked invalid.
   */
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;
  /**
   * Explicitly store traversal times for getEvidence().
   * Refrain from reading the data from GaussianProcesses to maintain abstraction.
   */
  std::unordered_map<Configuration, long, ConfigHash> _traversalTimes;

  Random _rng;
  /**
   * Stochastic model used for predictions.
   */
  GaussianCluster _gaussianCluster;
  /**
   * Function to generate neighbours with weight of given vector.
   */
  GaussianModelTypes::NeighbourFunction _neighbourFun;
  const size_t _maxEvidence;
  /**
   * Acquisition function used to predict informational gain.
   */
  const AcquisitionFunctionOption _predAcqFunction;
  /**
   * Number of latin-hypercube-samples used to find a evidence with high predicted acquisition.
   */
  const size_t _predNumLHSamples;

  bool _firstTuningPhase;
  /**
   * Iteration of last added evidence or reset.
   */
  size_t _currentIteration;
  /**
   * Iteration numbers get scaled down by a multiple of max evidence to ensure not all kernels become 0.
   */
  const double _iterationScale;
  /**
   * Number of evidence provided in current tuning phase.
   */
  size_t _currentNumEvidence;
  /**
   * Lowest time of evidence in current tuning phase.
   */
  long _currentOptimalTime;
  /**
   * Configuration with lowest time in current tuning phase.
   */
  FeatureVector _currentOptimalConfig;

  /**
   * FullSearch used in the first tuning phase.
   */
  FullSearch _fullSearch;
};

bool BayesianClusterSearch::tune(bool currentInvalid) {
  if (currentInvalid) {
    _invalidConfigs.insert(_currentConfig);
  }

  // in the first tuning phase do a full search
  if (_firstTuningPhase) {
    if (_fullSearch.tune(currentInvalid)) {
      // continue with full-search
      _currentConfig = _fullSearch.getCurrentConfiguration();
      return true;
    } else {
      // continue with GaussianCluster
      _firstTuningPhase = false;
      _currentNumEvidence = 0;
    }
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = FeatureVector();
    return false;
  }

  // no more tunings steps
  if (_currentNumEvidence >= _maxEvidence) {
    // select best config of current tuning phase
    _currentConfig = _currentOptimalConfig;

    return false;
  }

  // try to sample a valid vector which is expected to yield a good acquisition
  for (size_t i = 0; i < maxAttempts; ++i) {
    if (_currentAcquisitions.empty()) {
      sampleAcquisitions(_predNumLHSamples, _predAcqFunction);
    }

    // test best vectors until empty
    while (not _currentAcquisitions.empty()) {
      auto best = _encoder.convertFromCluster(_currentAcquisitions.back());

      if (_invalidConfigs.find(best) == _invalidConfigs.end()) {
        // valid config found!
        _currentConfig = best;
        return true;
      } else {
        // invalid: dispose
        _currentAcquisitions.pop_back();
      }
    }

    // No valid configuration. This should rarely happen.
    AutoPasLog(debug, "Tuning could not generate a valid configuration.");
  }

  utils::ExceptionHandler::exception("BayesianClusterSearch: Failed to sample an valid FeatureVector");
  _currentConfig = FeatureVector();
  return false;
}

void BayesianClusterSearch::sampleAcquisitions(size_t n, AcquisitionFunctionOption af) {
  // create n lhs samples
  auto continuousSamples = _encoder.lhsSampleFeatureCluster(n, _rng, _currentIteration * _iterationScale);

  // calculate all acquisitions
  _currentAcquisitions = _gaussianCluster.sampleOrderedByAcquisition(af, _neighbourFun, continuousSamples);
}

bool BayesianClusterSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerTraversalEstimatorOptions.size() == 1 and
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 1) and _dataLayoutOptions.size() == 1 and
         _newton3Options.size() == 1;
}

bool BayesianClusterSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerTraversalEstimatorOptions.empty() or
         (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or _dataLayoutOptions.empty() or
         _newton3Options.empty();
}

void BayesianClusterSearch::removeN3Option(Newton3Option badNewton3Option) {
  _fullSearch.removeN3Option(badNewton3Option);

  _newton3Options.erase(std::remove(_newton3Options.begin(), _newton3Options.end(), badNewton3Option),
                        _newton3Options.end());

  updateOptions();
  _currentAcquisitions.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

void BayesianClusterSearch::updateOptions() {
  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options,
                             *_cellSizeFactors, *_verletRebuildFrequencies);

  auto newRestrictions = _encoder.getDiscreteRestrictions();
  _gaussianCluster.setDimensions(std::vector<int>(newRestrictions.begin(), newRestrictions.end()));
}

}  // namespace autopas
