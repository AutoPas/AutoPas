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

#include "GaussianModel/GaussianCluster.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/FeatureVector.h"
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
   */
  constexpr static size_t continuousDims = 1;
  /**
   * Fixed noise
   */
  constexpr static double sigma = 0.01;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
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
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      size_t predNumLHSamples = 50, unsigned long seed = std::random_device()())
      : _containerOptionsSet(allowedContainerOptions),
        _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
        _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _currentConfig(),
        _invalidConfigs(),
        _rng(seed),
        _gaussianCluster(
            {static_cast<int>(allowedTraversalOptions.size()), static_cast<int>(allowedDataLayoutOptions.size()),
             static_cast<int>(allowedNewton3Options.size())},
            continuousDims, GaussianCluster::DistanceFunction::evidenceMatchingPDF, sigma, _rng),
        _neighbourFun([this](const Eigen::VectorXi &target) -> std::vector<Eigen::VectorXi> {
          return FeatureVector::neighboursManhattan1(target, _gaussianCluster.getDimensions());
        }),
        _maxEvidence(maxEvidence),
        _predAcqFunction(predAcqFunction),
        _predNumLHSamples(predNumLHSamples) {
    if (predNumLHSamples <= 0) {
      utils::ExceptionHandler::exception(
          "BayesianSearch: Number of samples used for predictions must be greater than 0!");
    }

    // generate all compatible container-traversal pairs
    for (auto traversal : allowedTraversalOptions) {
      auto compatibleContainerOptions = compatibleTraversals::allCompatibleContainers(traversal);
      std::set<ContainerOption> allowedAndCompatible;
      std::set_intersection(allowedContainerOptions.begin(), allowedContainerOptions.end(),
                            compatibleContainerOptions.begin(), compatibleContainerOptions.end(),
                            std::inserter(allowedAndCompatible, allowedAndCompatible.begin()));

      for (auto container : allowedAndCompatible) {
        _containerTraversalOptions.emplace_back(container, traversal);
      }
    }

    if (searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("BayesianClusterSearch: No valid configurations could be created.");
    }

    tune();
  }

  inline const Configuration &getCurrentConfiguration() const override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t) override {
    // encoded vector
    auto vec = _currentConfig.convertToCluster(_containerTraversalOptions, _dataLayoutOptions, _newton3Options);
    // time is converted to seconds, to big values may lead to errors in GaussianProcess. Time is also negated to
    // represent a maximization problem
    _gaussianCluster.addEvidence(vec, -time * secondsPerMicroseconds);
    _currentAcquisitions.clear();
  }

  inline void reset(size_t) override {
    _gaussianCluster.clear();
    _currentAcquisitions.clear();
    tune();
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
  std::vector<std::pair<ContainerOption, TraversalOption>> _containerTraversalOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  FeatureVector _currentConfig;
  /**
   * Currently sampled vectors and corresponding acquisition values.
   */
  std::vector<GaussianCluster::VectorPairDiscreteContinuous> _currentAcquisitions;
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
   * Function to generate neighbours of given vector
   */
  std::function<std::vector<Eigen::VectorXi>(const Eigen::VectorXi&)> _neighbourFun;
  const size_t _maxEvidence;
  /**
   * Acquisition function used to predict informational gain.
   */
  const AcquisitionFunctionOption _predAcqFunction;
  /**
   * Number of latin-hypercube-samples used to find a evidence with high predicted acquisition.
   */
  const size_t _predNumLHSamples;
};

bool BayesianClusterSearch::tune(bool currentInvalid) {
  if (currentInvalid) {
    _invalidConfigs.insert(_currentConfig);
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = FeatureVector();
    return false;
  }

  // no more tunings steps
  if (_gaussianCluster.numEvidence() >= _maxEvidence) {
    // select best config
    _currentConfig = FeatureVector::convertFromCluster(_gaussianCluster.getEvidenceMax(), _containerTraversalOptions,
                                                       _dataLayoutOptions, _newton3Options);
    AutoPasLog(debug, "Selected Configuration {}", _currentConfig.toString());

    // print graph of distances
    size_t sampleSize = _cellSizeFactors->isFinite() ? _cellSizeFactors->size() : _predNumLHSamples;
    auto continuousSamples = FeatureVector::lhsSampleFeatureContinuous(sampleSize, _rng, *_cellSizeFactors);
    auto vecToStringFun = [this] (const GaussianCluster::VectorPairDiscreteContinuous& vec) -> std::string {
      return FeatureVector::convertFromCluster(vec, _containerTraversalOptions, _dataLayoutOptions, _newton3Options).toString();
    };
    _gaussianCluster.logDebugGraph(_neighbourFun, continuousSamples, vecToStringFun);
    return false;
  }

  // try to sample a valid vector which is expected to yield a good acquisition
  for (size_t i = 0; i < maxAttempts; ++i) {
    if (_currentAcquisitions.empty()) {
      sampleAcquisitions(_predNumLHSamples, _predAcqFunction);
    }

    // test best vectors until empty
    while (not _currentAcquisitions.empty()) {
      auto best = FeatureVector::convertFromCluster(_currentAcquisitions.back(), _containerTraversalOptions,
                                                    _dataLayoutOptions, _newton3Options);

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
  auto continuousSamples = FeatureVector::lhsSampleFeatureContinuous(n, _rng, *_cellSizeFactors);

  // calculate all acquisitions
  _currentAcquisitions = _gaussianCluster.sampleOrderedByAcquisition(af, _neighbourFun, continuousSamples);
}

bool BayesianClusterSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerTraversalOptions.size() == 1 and (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 1) and
         _dataLayoutOptions.size() == 1 and _newton3Options.size() == 1;
}

bool BayesianClusterSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerTraversalOptions.empty() or (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or
         _dataLayoutOptions.empty() or _newton3Options.empty();
}

void BayesianClusterSearch::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(std::remove(_newton3Options.begin(), _newton3Options.end(), badNewton3Option),
                        _newton3Options.end());

  _gaussianCluster.setDimension(discreteNewtonDim, _newton3Options.size());
  _currentAcquisitions.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
