/**
 * @file SetSearchSpaceBasedTuningStrategy.cpp
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "SetSearchSpaceBasedTuningStrategy.h"

autopas::SetSearchSpaceBasedTuningStrategy::SetSearchSpaceBasedTuningStrategy(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
    const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options)
    : _allowedContainerOptions(allowedContainerOptions) {
  // sets search space and current config
  populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                      allowedLoadEstimatorOptions, allowedDataLayoutOptions, allowedNewton3Options);
}

autopas::SetSearchSpaceBasedTuningStrategy::SetSearchSpaceBasedTuningStrategy(
    std::set<Configuration> allowedConfigurations)
    : _allowedContainerOptions{}, _searchSpace(std::move(allowedConfigurations)) {
  for (const auto &config : _searchSpace) {
    _allowedContainerOptions.insert(config.container);
  }
  _currentConfig = _searchSpace.begin();
}

void autopas::SetSearchSpaceBasedTuningStrategy::populateSearchSpace(
    const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
    const std::set<TraversalOption> &allowedTraversalOptions,
    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
    const std::set<DataLayoutOption> &allowedDataLayoutOptions, const std::set<Newton3Option> &allowedNewton3Options) {
  // generate all potential configs
  for (const auto &containerOption : allowedContainerOptions) {
    // get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (const auto &cellSizeFactor : allowedCellSizeFactors)
      for (const auto &traversalOption : allowedAndApplicable) {
        // if load estimators are not applicable LoadEstimatorOption::none is returned.
        const std::set<LoadEstimatorOption> allowedAndApplicableLoadEstimators =
            loadEstimators::getApplicableLoadEstimators(containerOption, traversalOption, allowedLoadEstimatorOptions);
        for (const auto &loadEstimatorOption : allowedAndApplicableLoadEstimators) {
          for (const auto &dataLayoutOption : allowedDataLayoutOptions) {
            for (const auto &newton3Option : allowedNewton3Options) {
              _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, loadEstimatorOption,
                                   dataLayoutOption, newton3Option);
            }
          }
        }
      }
  }

  AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());

  if (_searchSpace.empty()) {
    autopas::utils::ExceptionHandler::exception("No valid configurations could be created.");
  }

  _currentConfig = _searchSpace.begin();
}

std::unordered_map<autopas::Configuration, long, autopas::ConfigHash>::const_iterator
autopas::SetSearchSpaceBasedTuningStrategy::getOptimum(
    const std::unordered_map<Configuration, long, ConfigHash> &searchSpace) {
  // find mapping with the smallest second (=runtime)
  return std::min_element(
      searchSpace.begin(), searchSpace.end(),
      [](std::pair<Configuration, long> a, std::pair<Configuration, long> b) -> bool { return a.second < b.second; });
}

void autopas::SetSearchSpaceBasedTuningStrategy::removeN3Option(Newton3Option badNewton3Option) {
  for (auto ssIter = _searchSpace.begin(); ssIter != _searchSpace.end();) {
    if (ssIter->newton3 == badNewton3Option) {
      // change current config to the next non-deleted
      if (ssIter == _currentConfig) {
        ssIter = _searchSpace.erase(ssIter);
        _currentConfig = ssIter;
      } else {
        ssIter = _searchSpace.erase(ssIter);
      }
    } else {
      ++ssIter;
    }
  }

  if (this->searchSpaceIsEmpty()) {
    utils::ExceptionHandler::exception(
        "Removing all configurations with Newton 3 {} caused the search space to be empty!", badNewton3Option);
  }
}

std::set<autopas::ContainerOption> autopas::SetSearchSpaceBasedTuningStrategy::getAllowedContainerOptions() const {
  return _allowedContainerOptions;
}

const autopas::Configuration &autopas::SetSearchSpaceBasedTuningStrategy::getCurrentConfiguration() const {
  return *_currentConfig;
}

bool autopas::SetSearchSpaceBasedTuningStrategy::searchSpaceIsTrivial() const { return _searchSpace.size() == 1; }

bool autopas::SetSearchSpaceBasedTuningStrategy::searchSpaceIsEmpty() const { return _searchSpace.empty(); }

bool autopas::SetSearchSpaceBasedTuningStrategy::smoothedHomogeneityAndMaxDensityNeeded() const { return false; }
