/**
 * @file RandomSearch.h
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "RandomSearch.h"

autopas::RandomSearch::RandomSearch(const std::set<ContainerOption> &allowedContainerOptions,
                                    const autopas::NumberSet<double> &allowedCellSizeFactors,
                                    const std::set<TraversalOption> &allowedTraversalOptions,
                                    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                    const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                    const std::set<Newton3Option> &allowedNewton3Options, size_t maxEvidence,
                                    unsigned long seed)
    : _containerOptions(allowedContainerOptions),
      _traversalOptions(allowedTraversalOptions),
      _loadEstimatorOptions(allowedLoadEstimatorOptions),
      _dataLayoutOptions(allowedDataLayoutOptions),
      _newton3Options(allowedNewton3Options),
      _cellSizeFactors(allowedCellSizeFactors.clone()),
      _currentConfig(),
      _invalidConfigs(),
      _rng(seed),
      _maxEvidence(maxEvidence),
      _searchSpaceSizeNoCSF(0),
      _numTestedConfigsNoCSF(0) {
  if (searchSpaceIsEmpty()) {
    autopas::utils::ExceptionHandler::exception("RandomSearch: No valid configurations could be created.");
  }

  // Allow invalid container-traversal-combinations, because so does tune().
  for (auto container : _containerOptions) {
    for (auto traversal : _traversalOptions) {
      auto applicableLoadEstimators =
          loadEstimators::getApplicableLoadEstimators(container, traversal, _loadEstimatorOptions);
      _searchSpaceSizeNoCSF += applicableLoadEstimators.size();
    }
  }
  _searchSpaceSizeNoCSF *= _dataLayoutOptions.size() * _newton3Options.size();

  // if fewer than _maxEvidence configurations exist in total, search through all of them and then stop.
  _maxEvidence = std::min(_searchSpaceSizeNoCSF, _maxEvidence);
  tune();
}

bool autopas::RandomSearch::tune(bool currentInvalid) {
  // We can overwrite the cellSizeFactor, because _currentConfig will be changed anyways.
  if (currentInvalid) {
    _currentConfig.cellSizeFactor = 1;
    if (_invalidConfigs.find(_currentConfig) == _invalidConfigs.end()) {
      ++_numTestedConfigsNoCSF;
    }
    _invalidConfigs.emplace(_currentConfig);
  }

  if (searchSpaceIsEmpty()) {
    // no valid configuration
    _currentConfig = Configuration();
    return false;
  }

  if (_traversalTimes.size() >= _maxEvidence or _numTestedConfigsNoCSF >= _searchSpaceSizeNoCSF) {
    // select best config
    selectOptimalConfiguration();
    return false;
  }

  // needed to test if a new configuration has been found for infinite cellSizeFactors.
  Configuration testingConfig;
  do {
    // select random config
    _currentConfig.container = _rng.pickRandom(_containerOptions);
    _currentConfig.cellSizeFactor = _cellSizeFactors->getRandom(_rng);
    _currentConfig.traversal = _rng.pickRandom(_traversalOptions);
    auto applicableLoadEstimators = loadEstimators::getApplicableLoadEstimators(
        _currentConfig.container, _currentConfig.traversal, _loadEstimatorOptions);
    _currentConfig.loadEstimator = _rng.pickRandom(applicableLoadEstimators);
    _currentConfig.dataLayout = _rng.pickRandom(_dataLayoutOptions);
    _currentConfig.newton3 = _rng.pickRandom(_newton3Options);

    testingConfig = Configuration(_currentConfig);
    if (not _cellSizeFactors->isFinite()) {
      testingConfig.cellSizeFactor = 1;
    }
  } while (_invalidConfigs.find(testingConfig) != _invalidConfigs.end() or
           _traversalTimes.find(_currentConfig) != _traversalTimes.end());
  return true;
}

void autopas::RandomSearch::selectOptimalConfiguration() {
  if (searchSpaceIsTrivial()) {
    _currentConfig.container = *_containerOptions.begin();
    _currentConfig.cellSizeFactor = _cellSizeFactors->getMin();
    _currentConfig.traversal = *_traversalOptions.begin();
    auto applicableLoadEstimators = loadEstimators::getApplicableLoadEstimators(
        _currentConfig.container, _currentConfig.traversal, _loadEstimatorOptions);
    _currentConfig.loadEstimator = *applicableLoadEstimators.begin();
    _currentConfig.dataLayout = *_dataLayoutOptions.begin();
    _currentConfig.newton3 = *_newton3Options.begin();
    AutoPasLog(debug, "Selected optimal configuration {}", _currentConfig.toString());
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "RandomSearch: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  auto optimum = std::min_element(_traversalTimes.begin(), _traversalTimes.end(),
                                  [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                                    return a.second < b.second;
                                  });

  _currentConfig = optimum->first;

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig.toString());
}

bool autopas::RandomSearch::searchSpaceIsTrivial() const {
  if (searchSpaceIsEmpty()) {
    return false;
  }

  // if no load estimators are specified, none is automatically added.
  return _containerOptions.size() == 1 and (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 1) and
         _traversalOptions.size() == 1 and _loadEstimatorOptions.size() <= 1 and _dataLayoutOptions.size() == 1 and
         _newton3Options.size() == 1;
}

bool autopas::RandomSearch::searchSpaceIsEmpty() const {
  // if one enum is empty return true.
  return _containerOptions.empty() or (_cellSizeFactors->isFinite() and _cellSizeFactors->size() == 0) or
         _traversalOptions.empty() or _dataLayoutOptions.empty() or _newton3Options.empty();
}

void autopas::RandomSearch::removeN3Option(Newton3Option badNewton3Option) {
  _newton3Options.erase(badNewton3Option);

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

const autopas::Configuration &autopas::RandomSearch::getCurrentConfiguration() const { return _currentConfig; }

void autopas::RandomSearch::addEvidence(long time, size_t iteration) {
  // Count the number of valid, tested configurations to test if there are still untested configurations that could be
  // valid.
  // Ignore CellSizeFactors because infinite CSFs make this test difficult.
  auto testingConfig = _currentConfig;
  testingConfig.cellSizeFactor = 1;
  if (_traversalTimes.find(testingConfig) == _traversalTimes.end()) {
    ++_numTestedConfigsNoCSF;
  }

  _traversalTimes[_currentConfig] = time;
}

long autopas::RandomSearch::getEvidence(autopas::Configuration configuration) const {
  return static_cast<long>(_traversalTimes.at(configuration));
}

void autopas::RandomSearch::reset(size_t iteration) {
  _traversalTimes.clear();
  _invalidConfigs.clear();
  _numTestedConfigsNoCSF = 0;
  tune();
}

std::set<autopas::ContainerOption> autopas::RandomSearch::getAllowedContainerOptions() const {
  return _containerOptions;
}
