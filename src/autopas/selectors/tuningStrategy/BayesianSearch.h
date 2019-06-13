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
#include "autopas/utils/DoubleSet.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given
 * feature next.
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
   * @param expAcqFunction acquisition function used for exploration.
   * @param expNumSamples number of samples used for exploration.
   */
  BayesianSearch(const std::set<ContainerOption>& allowedContainerOptions = allContainerOptions,
                 const std::set<TraversalOption>& allowedTraversalOptions = allTraversalOptions,
                 const std::set<DataLayoutOption>& allowedDataLayoutOptions = allDataLayoutOptions,
                 const std::set<Newton3Option>& allowedNewton3Options = allNewton3Options,
                 const DoubleSet& allowedCellSizeFactors = DoubleInterval(0., 2.),
                 AcquisitionFunction expAcqFunction = lcb, size_t expNumSamples = 1000)
      : _enumOptions(4),
        _containerOptions(allowedContainerOptions),
        _traversalOptions(allowedTraversalOptions),
        _dataLayoutOptions(allowedDataLayoutOptions),
        _newton3Options(allowedNewton3Options),
        _cellSizeFactors(allowedCellSizeFactors.clone()),
        _currentConfig(),
        _gp(0., std::vector<double>(_enumOptions.size() + 1, 1.), 0.001),
        _expAcqFunction(expAcqFunction),
        _expNumSamples(expNumSamples),
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
  void updateEnumOptions() {
    _enumOptions = {EnumFeature::set2Vector(_containerOptions), EnumFeature::set2Vector(_traversalOptions),
                    EnumFeature::set2Vector(_dataLayoutOptions), EnumFeature::set2Vector(_newton3Options)};
  }

  /**
   * Set the current configuration with
   * information from a FeatureVector
   */
  void setConfig(FeatureVector feature);
  /**
   * Get the current configuration as
   * a FeatureVector.
   */
  FeatureVector configAsFeature();

  std::vector<std::vector<EnumFeature>> _enumOptions;
  std::set<ContainerOption> _containerOptions;
  std::set<TraversalOption> _traversalOptions;
  std::set<DataLayoutOption> _dataLayoutOptions;
  std::set<Newton3Option> _newton3Options;
  std::unique_ptr<DoubleSet> _cellSizeFactors;

  Configuration _currentConfig;
  GaussianProcess _gp;
  AcquisitionFunction _expAcqFunction;
  size_t _expNumSamples;
  std::default_random_engine _rng;
};

bool BayesianSearch::tune() {
  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = Configuration();
    return false;
  }

  std::vector<FeatureVector> samples(_expNumSamples);

  // sample from all enums
  for (auto& enumOption : _enumOptions) {
    FeatureVector::lhsAddFeature(samples, enumOption, _rng);
  }
  // sample from cellSizeFactors
  FeatureVector::lhsAddFeature(samples, *_cellSizeFactors, _rng);

  // sample minimum of acquisition function
  FeatureVector best = _gp.sampleAquisitionMin(_expAcqFunction, samples);
  setConfig(best);

  return true;
}

bool BayesianSearch::searchSpaceIsTrivial() {
  bool multiOptions = false;

  // check all enums
  for (auto& _enumOption : _enumOptions) {
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
  for (auto& _enumOption : _enumOptions) {
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
  ContainerOption co = static_cast<ContainerOption>(dynamic_cast<EnumFeature&>(feature[0]).getValue());
  TraversalOption to = static_cast<TraversalOption>(dynamic_cast<EnumFeature&>(feature[1]).getValue());
  DataLayoutOption dlo = static_cast<DataLayoutOption>(dynamic_cast<EnumFeature&>(feature[2]).getValue());
  Newton3Option n3o = static_cast<Newton3Option>(dynamic_cast<EnumFeature&>(feature[3]).getValue());
  double csf = dynamic_cast<DoubleFeature&>(feature[4]).getValue();

  _currentConfig = Configuration(co, to, dlo, n3o, csf);
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
