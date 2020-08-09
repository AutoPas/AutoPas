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
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LoadEstimators.h"
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
   * Dimension corresponding to the newton3 option in the discrete tuple.
   */
  constexpr static size_t discreteNewtonDim = 2;
  /**
   * Dimensions of the continuous tuples.
   * FeatureVector continuous dimensions + 1 (time)
   */
  constexpr static size_t continuousDims = FeatureVector::featureSpaceContinuousDims + 1;
  /**
   * Fixed noise
   */
  constexpr static double sigma = 0.01;

  /**
   * Iterations get scaled down by a multiple of max evidence to ensure not all kernels become 0.
   */
  constexpr static double iterationScalePerMaxEvidence = 5.;

  /**
   * Distance between evidence is expected not to exceed a upper threshold.
   * We chose the threshold 7 because: exp(-r) < 0.1% for r > 7 but this
   * can be increased if necessary.
   */
  constexpr static double suggestedMaxDistance = 7.;

  /**
   * Number of cellSizeFactors sampled for FullSearch.
   */
  constexpr static size_t cellSizeFactorSampleSize = 1;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumLHSamples number of latin-hypercube-samples used to find a evidence with high predicted acquisition
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  explicit BayesianClusterSearch(
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      size_t predNumLHSamples = 50, unsigned long seed = std::random_device()())
      : _containerOptionsSet(allowedContainerOptions),
        _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
        _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _encoder(),
        _currentConfig(),
        _invalidConfigs(),
        _rng(seed),
        _gaussianCluster({}, continuousDims, GaussianCluster::WeightFunction::evidenceMatchingScaledProbabilityGM,
                         sigma, _rng),
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
        _fullSearch(allowedContainerOptions, allowedCellSizeFactors.uniformSampleSet(cellSizeFactorSampleSize, _rng),
                    allowedTraversalOptions, allowedLoadEstimatorOptions, allowedDataLayoutOptions,
                    allowedNewton3Options) {
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

    _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options);
    _gaussianCluster.setDimensions({static_cast<int>(_containerTraversalEstimatorOptions.size()),
                                    static_cast<int>(_dataLayoutOptions.size()),
                                    static_cast<int>(_newton3Options.size())});
    _gaussianCluster.setVectorToStringFun(
        [this](const GaussianModelTypes::VectorPairDiscreteContinuous &vec) -> std::string {
          return _encoder.convertFromClusterWithIteration(vec).toString();
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
    auto vec = _encoder.convertToClusterWithIteration(_currentConfig, iteration * _iterationScale);
    // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
    // represent a maximization problem
    _gaussianCluster.addEvidence(vec, -time * secondsPerMicroseconds);

    _currentIteration = iteration;
    ++_currentNumEvidence;
    _currentAcquisitions.clear();
  }

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
    } else {
      tune();
    }
  }

  inline bool tune(bool currentInvalid = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptionsSet; }

  inline bool searchSpaceIsTrivial() const override;

  inline bool searchSpaceIsEmpty() const override;

 private:
  /**
   * Generate n samples and predict their corresponding
   * acquisition function. The result is sorted depending on
   * the given acquisition function.
   * @param n numSamples
   * @param af acquisition function
   */
  inline void sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
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
   * We scale current iteration down such that kernels do not become 0.
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
    AutoPasLog(debug, "Selected Configuration {}", _currentConfig.toString());

    return false;
  }

  // try to sample a valid vector which is expected to yield a good acquisition
  for (size_t i = 0; i < maxAttempts; ++i) {
    if (_currentAcquisitions.empty()) {
      sampleAcquisitions(_predNumLHSamples, _predAcqFunction);
    }

    // test best vectors until empty
    while (not _currentAcquisitions.empty()) {
      auto best = _encoder.convertFromClusterWithIteration(_currentAcquisitions.back());

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

  utils::ExceptionHandler::exception("BayesianSearch: Failed to sample an valid FeatureVector");
  _currentConfig = FeatureVector();
  return false;
}

void BayesianClusterSearch::sampleAcquisitions(size_t n, AcquisitionFunctionOption af) {
  // create n lhs samples
  auto continuousSamples = FeatureVector::lhsSampleFeatureContinuousWithIteration(n, _rng, *_cellSizeFactors,
                                                                                  _currentIteration * _iterationScale);

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
  _encoder.setAllowedOptions(_containerTraversalEstimatorOptions, _dataLayoutOptions, _newton3Options);

  _gaussianCluster.setDimensions({static_cast<int>(_containerTraversalEstimatorOptions.size()),
                                  static_cast<int>(_dataLayoutOptions.size()),
                                  static_cast<int>(_newton3Options.size())});
  _currentAcquisitions.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}
}  // namespace autopas
