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
   *
   */
  constexpr static double secondsPerMicroseconds = 1. / 1000000.;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumLHSamples number of samples used for prediction while tuning.
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  BayesianClusterSearch(const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
                        const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
                        const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
                        const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
                        const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(),
                        size_t maxEvidence = 10,
                        AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::lowerConfidenceBound,
                        size_t predNumLHSamples = 50, unsigned long seed = std::random_device()())
      : _containerOptionsSet(allowedContainerOptions),
        _containerOptions(allowedContainerOptions.begin(), allowedContainerOptions.end()),
        _traversalOptions(allowedTraversalOptions.begin(), allowedTraversalOptions.end()),
        _dataLayoutOptions(allowedDataLayoutOptions.begin(), allowedDataLayoutOptions.end()),
        _newton3Options(allowedNewton3Options.begin(), allowedNewton3Options.end()),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _traversalContainerMap(),
        _currentConfig(),
        _invalidConfigs(),
        _rng(seed),
        _dimRestrictions({static_cast<int>(allowedTraversalOptions.size()),
                          static_cast<int>(allowedDataLayoutOptions.size()),
                          static_cast<int>(allowedNewton3Options.size())}),
        _gaussianCluster(_dimRestrictions, 1, 0.01, _rng),
        _maxEvidence(maxEvidence),
        _predAcqFunction(predAcqFunction),
        _predNumLHSamples(predNumLHSamples) {
    if (predNumLHSamples <= 0) {
      utils::ExceptionHandler::exception(
          "BayesianSearch: Number of samples used for predictions must be greater than 0!");
    }

    // map each travesal to a container
    for (auto it = _traversalOptions.begin(); it != _traversalOptions.end();) {
      auto compatibleContainerOptions = compatibleTraversals::allCompatibleContainers(*it);
      std::set<ContainerOption> allowedAndCompatible;
      std::set_intersection(allowedContainerOptions.begin(), allowedContainerOptions.end(),
                            compatibleContainerOptions.begin(), compatibleContainerOptions.end(),
                            std::inserter(allowedAndCompatible, allowedAndCompatible.begin()));
      if (allowedAndCompatible.empty()) {
        // no compatible container found
        it = _traversalOptions.erase(it);
      } else {
        _traversalContainerMap[*it] = *allowedAndCompatible.begin();
        ++it;
      }
    }

    tune();
  }

  inline const Configuration &getCurrentConfiguration() const override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time) override {
    // encoded vector
    auto vec = _currentConfig.clusterEncode(_traversalOptions, _dataLayoutOptions, _newton3Options);
    // time is converted to seconds, to big values may lead to errors in GaussianProcess
    _gaussianCluster.addEvidence(vec.first, vec.second, time * secondsPerMicroseconds);
    _currentAcquisitions.clear();
  }

  inline void reset() override {
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
   * given acquistion function.
   * @param n numSamples
   * @param af acquisition function
   */
  inline void sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<ContainerOption> _containerOptions;
  std::vector<TraversalOption> _traversalOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  std::map<TraversalOption, ContainerOption> _traversalContainerMap;

  FeatureVector _currentConfig;
  std::vector<GaussianCluster::VectorAcquisition> _currentAcquisitions;
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;

  Random _rng;
  std::vector<int> _dimRestrictions;
  GaussianCluster _gaussianCluster;
  size_t _maxEvidence;
  AcquisitionFunctionOption _predAcqFunction;
  size_t _predNumLHSamples;
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
    _currentConfig = FeatureVector::clusterDecode(_gaussianCluster.getEvidenceMin(), _traversalOptions,
                                                  _dataLayoutOptions, _newton3Options);
    _currentConfig.container = _traversalContainerMap[_currentConfig.traversal];
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
      auto best = FeatureVector::clusterDecode(_currentAcquisitions[0].first, _traversalOptions, _dataLayoutOptions,
                                               _newton3Options);
      best.container = _traversalContainerMap[best.traversal];

      if (_invalidConfigs.find(best) == _invalidConfigs.end()) {
        // valid config found!
        _currentConfig = best;
        return true;
      } else {
        // invalid: dispose
        _currentAcquisitions.erase(_currentAcquisitions.begin());
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

  auto neighbourFun = [this](Eigen::VectorXi target) -> std::vector<Eigen::VectorXi> {
    return FeatureVector::neighboursManhattan1(target, _dimRestrictions);
  };

  // calculate all acquisitions
  _currentAcquisitions = _gaussianCluster.sampleAcquisition(af, neighbourFun, continuousSamples);

  // sort vectors by acquisition
  switch (af) {
    case AcquisitionFunctionOption::mean:
    case AcquisitionFunctionOption::lowerConfidenceBound:
    case AcquisitionFunctionOption::upperConfidenceBound:
      // min estimated value first
      std::sort(_currentAcquisitions.begin(), _currentAcquisitions.end(),
                [](const GaussianCluster::VectorAcquisition &va1, const GaussianCluster::VectorAcquisition &va2) {
                  return va1.second < va2.second;
                });
      break;
    case AcquisitionFunctionOption::variance:
    case AcquisitionFunctionOption::expectedDecrease:
    case AcquisitionFunctionOption::probabilityOfDecrease:
      // max acquisition first
      std::sort(_currentAcquisitions.begin(), _currentAcquisitions.end(),
                [](const GaussianCluster::VectorAcquisition &va1, const GaussianCluster::VectorAcquisition &va2) {
                  return va1.second > va2.second;
                });
      break;
  }
}

bool BayesianClusterSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  return _containerOptions.size() == 1 and (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 1) and
         _traversalOptions.size() == 1 and _dataLayoutOptions.size() == 1 and _newton3Options.size() == 1;
}

bool BayesianClusterSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true
  return _containerOptions.empty() or (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 0) or
         _traversalOptions.empty() or _dataLayoutOptions.empty() or _newton3Options.empty();
}

void BayesianClusterSearch::removeN3Option(Newton3Option badNewton3Option) {
  std::remove(_newton3Options.begin(), _newton3Options.end(), badNewton3Option);
  _currentAcquisitions.clear();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

}  // namespace autopas
