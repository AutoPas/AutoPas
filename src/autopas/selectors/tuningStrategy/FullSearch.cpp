/**
 * @file FullSearch.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "FullSearch.h"

autopas::FullSearch::FullSearch(const std::set<ContainerOption> &allowedContainerOptions,
                                const std::set<double> &allowedCellSizeFactors,
                                const std::set<TraversalOption> &allowedTraversalOptions,
                                const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                const std::set<Newton3Option> &allowedNewton3Options)
    : SetSearchSpaceBasedTuningStrategy(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                                        allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options) {}

autopas::FullSearch::FullSearch(std::set<Configuration> allowedConfigurations)
    : SetSearchSpaceBasedTuningStrategy(std::move(allowedConfigurations)) {}

bool autopas::FullSearch::tune(bool) {
  // repeat as long as traversals are not applicable or we run out of configs
  ++_currentConfig;
  if (_currentConfig == _searchSpace.end()) {
    selectOptimalConfiguration();
    return false;
  }

  return true;
}

void autopas::FullSearch::selectOptimalConfiguration() {
  if (_searchSpace.size() == 1) {
    _currentConfig = _searchSpace.begin();
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "FullSearch: Trying to determine fastest configuration without any measurements! "
        "Either selectOptimalConfiguration was called too early or no applicable configurations were found");
  }

  const auto optimum = std::min_element(_traversalTimes.begin(), _traversalTimes.end(),
                                        [](std::pair<Configuration, size_t> a,
                                           std::pair<Configuration, size_t> b) -> bool { return a.second < b.second; });

  _currentConfig = _searchSpace.find(optimum->first);
  // sanity check
  if (_currentConfig == _searchSpace.end()) {
    autopas::utils::ExceptionHandler::exception(
        "FullSearch: Optimal configuration not found in list of configurations!");
  }
}

void autopas::FullSearch::addEvidence(long time, size_t iteration) { _traversalTimes[*_currentConfig] = time; }

long autopas::FullSearch::getEvidence(autopas::Configuration configuration) const {
  return _traversalTimes.at(configuration);
}

const autopas::Configuration &autopas::FullSearch::getCurrentConfiguration() const { return *_currentConfig; }

void autopas::FullSearch::reset(size_t iteration) {
  _traversalTimes.clear();
  _currentConfig = _searchSpace.begin();
}
