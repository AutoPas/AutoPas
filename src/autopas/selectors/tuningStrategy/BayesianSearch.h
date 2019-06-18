/**
 * @file BayesianSearch.h
 * @author Jan Nguyen
 * @date 12.06.19
 */

#pragma once

#include <set>
#include <sstream>
#include "GaussianProcess.h"
#include "TuningStrategyInterface.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given
 * feature next. @todo find better conversion to featureVector
 */
class BayesianSearch : public TuningStrategyInterface {
 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumSamples number of samples used for prediction while tuning.
   * @param maxEvidences stop tuning after given number evidences provided.
   * @param lastAcqFunction acquisition function used for prediction of last tuning step.
   * @param lastNumSamples number of samples used for prediction of last tuning step.
   */
  BayesianSearch(const std::set<ContainerOption> &allowedContainerOptions = allContainerOptions,
                 const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(0., 2.),
                 const std::set<TraversalOption> &allowedTraversalOptions = allTraversalOptions,
                 const std::set<DataLayoutOption> &allowedDataLayoutOptions = allDataLayoutOptions,
                 const std::set<Newton3Option> &allowedNewton3Options = allNewton3Options,
                 AcquisitionFunction predAcqFunction = lcb, size_t predNumSamples = 1000, size_t maxEvidences = 10,
                 AcquisitionFunction lastAcqFunction = ucb, size_t lastNumSamples = 100000)
      : _enumOptions(4),
        _containerOptions(allowedContainerOptions),
        _traversalOptions(allowedTraversalOptions),
        _dataLayoutOptions(allowedDataLayoutOptions),
        _newton3Options(allowedNewton3Options),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _currentConfig(),
        _gp(0., std::vector<double>(_enumOptions.size() + 1, 1.), 0.001),
        _predAcqFunction(predAcqFunction),
        _predNumSamples(predNumSamples),
        _maxEvidences(maxEvidences),
        _lastAcqFunction(lastAcqFunction),
        _lastNumSamples(lastNumSamples),
        _rng(std::random_device()()) {
    /// @todo setting hyperparameters

    updateEnumOptions();
    tune();
  }

  inline Configuration getCurrentConfiguration() override { return _currentConfig; }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time) override { _gp.addEvidence(configAsFeature(), time); }

  inline void reset() override {
    _gp.clear();
    tune();
  }

  inline bool tune() override;

  inline std::set<ContainerOption> getAllowedContainerOptions() override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() override;

  inline bool searchSpaceIsEmpty() override;

 private:
  /**
   * Update search space for enums
   */
  inline void updateEnumOptions() {
    _enumOptions = {EnumFeature::set2Vector(_containerOptions), EnumFeature::set2Vector(_traversalOptions),
                    EnumFeature::set2Vector(_dataLayoutOptions), EnumFeature::set2Vector(_newton3Options)};
  }

  /**
   * Generate n samples and predict their corresponding
   * acquisition function. Return the FeatureVector which
   * generates the smallest value.
   * @param n numSamples
   * @param af acquisition function
   */
  inline FeatureVector sampleOptimalFeatureVector(size_t n, AcquisitionFunction af);

  /**
   * Set the current configuration with
   * information from a FeatureVector
   */
  inline void setConfig(FeatureVector feature);
  /**
   * Get the current configuration as
   * a FeatureVector.
   */
  inline FeatureVector configAsFeature();

  std::vector<std::vector<EnumFeature>> _enumOptions;
  std::set<ContainerOption> _containerOptions;
  std::set<TraversalOption> _traversalOptions;
  std::set<DataLayoutOption> _dataLayoutOptions;
  std::set<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  Configuration _currentConfig;
  GaussianProcess _gp;
  AcquisitionFunction _predAcqFunction;
  size_t _predNumSamples;
  size_t _maxEvidences;
  AcquisitionFunction _lastAcqFunction;
  size_t _lastNumSamples;
  std::default_random_engine _rng;
};

bool BayesianSearch::tune() {
  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = Configuration();
    return false;
  }

  if (_gp.numEvidences() >= _maxEvidences) {
    // select predicted best config
    setConfig(sampleOptimalFeatureVector(_lastNumSamples, _lastAcqFunction));
    return false;
  }

  // select predicted best config for tuning
  setConfig(sampleOptimalFeatureVector(_predNumSamples, _predAcqFunction));
  return true;
}

FeatureVector BayesianSearch::sampleOptimalFeatureVector(size_t n, AcquisitionFunction af) {
  std::vector<FeatureVector> samples;
  while (samples.empty()) {
    samples.resize(n);

    // sample from all enums
    for (auto &enumOption : _enumOptions) {
      FeatureVector::lhsAddFeature(samples, enumOption, _rng);
    }
    // sample from cellSizeFactors
    FeatureVector::lhsAddFeature(samples, *_cellSizeFactors, _rng);

    // remove all invalid samples
    for (auto it = samples.begin(); it != samples.end();) {
      auto allContainerTraversals = compatibleTraversals::allCompatibleTraversals(
          static_cast<ContainerOption>(dynamic_cast<EnumFeature &>((*it)[0]).getValue()));
      auto traversalOption = static_cast<TraversalOption>(dynamic_cast<EnumFeature &>((*it)[1]).getValue());
      if (allContainerTraversals.find(traversalOption) == allContainerTraversals.end()) {
        it = samples.erase(it);
      } else {
        ++it;
      }
    }
  }

  // sample minimum of acquisition function
  return _gp.sampleAquisitionMin(af, samples);
}

bool BayesianSearch::searchSpaceIsTrivial() {
  bool multiOptions = false;

  // check all enums
  for (auto &_enumOption : _enumOptions) {
    if (_enumOption.size() == 0) {
      // no option allowed
      return false;
    } else if (_enumOption.size() > 1) {
      // multiple options allowed
      multiOptions = true;
    }
  }

  // check csf
  if (_cellSizeFactors->isFinite()) {
    if (_cellSizeFactors->size() == 0) {
      // no option allowed
      return false;
    } else if (_cellSizeFactors->size() > 1) {
      // multiple options allowed
      multiOptions = true;
    }
  } else {
    // multiple options allowed
    multiOptions = true;
  }

  return not multiOptions;
}

bool BayesianSearch::searchSpaceIsEmpty() {
  // if one enum is empty return true
  for (auto &_enumOption : _enumOptions) {
    if (_enumOption.size() == 0) return true;
  }

  // if csf is empty return true
  if (_cellSizeFactors->isFinite() && _cellSizeFactors->size() == 0) {
    return true;
  }

  return false;
}

void BayesianSearch::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(badNewton3Option);
  updateEnumOptions();

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

void BayesianSearch::setConfig(FeatureVector feature) {
  ContainerOption co = static_cast<ContainerOption>(dynamic_cast<EnumFeature &>(feature[0]).getValue());
  TraversalOption to = static_cast<TraversalOption>(dynamic_cast<EnumFeature &>(feature[1]).getValue());
  DataLayoutOption dlo = static_cast<DataLayoutOption>(dynamic_cast<EnumFeature &>(feature[2]).getValue());
  Newton3Option n3o = static_cast<Newton3Option>(dynamic_cast<EnumFeature &>(feature[3]).getValue());
  double csf = dynamic_cast<DoubleFeature &>(feature[4]).getValue();

  _currentConfig = Configuration(co, csf, to, dlo, n3o);
}

FeatureVector BayesianSearch::configAsFeature() {
  FeatureVector result;
  result.addFeature(EnumFeature(_currentConfig.container));
  result.addFeature(EnumFeature(_currentConfig.traversal));
  result.addFeature(EnumFeature(_currentConfig.dataLayout));
  result.addFeature(EnumFeature(_currentConfig.newton3));
  result.addFeature(_currentConfig.cellSizeFactor);

  return result;
}

}  // namespace autopas
