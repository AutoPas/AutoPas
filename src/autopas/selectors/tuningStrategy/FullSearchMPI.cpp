/**
 * @file FullSearchMPI.cpp
 * @author W. Thieme
 * @date 29.05.2020
 */

#include "FullSearchMPI.h"

namespace autopas {

FullSearchMPI::FullSearchMPI(const std::set<ContainerOption> &allowedContainerOptions,
                             const std::set<double> &allowedCellSizeFactors,
                             const std::set<TraversalOption> &allowedTraversalOptions,
                             const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                             const std::set<Newton3Option> &allowedNewton3Options, AutoPas_MPI_Comm comm)
    : _containerOptions(allowedContainerOptions),
      _optimalConfig(Configuration()),
      _configurationCommunicator(AutoPasConfigurationCommunicator()),
      _localOptimalTime(0),
      _request(AUTOPAS_MPI_REQUEST_NULL),
      _allLocalConfigsTested(false),
      _allGlobalConfigsTested(false),
      _autopasMPICommunicator(comm) {
  // Count the total number of allowed and applicable configurations
  // Currently every rank does this independently and without communication, just based on rank and commSize
  int totalNumConfigs = 0;
  for (auto &containerOption : allowedContainerOptions) {
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));
    totalNumConfigs += allowedCellSizeFactors.size() * allowedAndApplicable.size() * allowedDataLayoutOptions.size() *
                       allowedNewton3Options.size();
  }

  int worldSize;
  AutoPas_MPI_Comm_size(_autopasMPICommunicator, &worldSize);

  int worldRank;
  AutoPas_MPI_Comm_rank(_autopasMPICommunicator, &worldRank);

  // divide search space into worldSize many blocks.
  const int blockSize = totalNumConfigs / worldSize;
  const int remainder = totalNumConfigs % worldSize;
  int start = blockSize * worldRank;
  int startNext = blockSize * (worldRank + 1);
  // distribute the remaining configurations over the lower ranks.
  if (worldRank < remainder) {
    start += worldRank;
    startNext += worldRank + 1;
  } else {
    start += remainder;
    startNext += remainder;
  }

  if (worldRank == 0) {
    AutoPasLog(debug, "number of ranks: {}", worldSize);
    AutoPasLog(debug, "total number of possible configurations: {}", totalNumConfigs);
  }
  populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                      allowedDataLayoutOptions, allowedNewton3Options, start, startNext);

  if (_searchSpace.empty()) {
    populateSearchSpace(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                        allowedDataLayoutOptions, allowedNewton3Options, 0, totalNumConfigs);
    if (_searchSpace.empty()) {
      utils::ExceptionHandler::exception("No valid configuration could be generated");
    }
  }
  AutoPasLog(debug, "Points in search space: {}", _searchSpace.size());
  _tuningConfig = _searchSpace.begin();
}

void FullSearchMPI::populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                        const std::set<double> &allowedCellSizeFactors,
                                        const std::set<TraversalOption> &allowedTraversalOptions,
                                        const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                        const std::set<Newton3Option> &allowedNewton3Options, const int start,
                                        const int startNext) {
  // Index to test which configurations apply to the current rank
  int i = 0;
  // Generate all potential configs
  for (auto &containerOption : allowedContainerOptions) {
    // Get all traversals of the container and restrict them to the allowed ones
    const std::set<TraversalOption> &allContainerTraversals =
        compatibleTraversals::allCompatibleTraversals(containerOption);
    std::set<TraversalOption> allowedAndApplicable;
    std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                          allContainerTraversals.begin(), allContainerTraversals.end(),
                          std::inserter(allowedAndApplicable, allowedAndApplicable.begin()));

    for (auto &cellSizeFactor : allowedCellSizeFactors) {
      for (auto &traversalOption : allowedAndApplicable) {
        for (auto &dataLayoutOption : allowedDataLayoutOptions) {
          for (auto &newton3Option : allowedNewton3Options) {
            if (i < start) {
              ++i;
              continue;
            }
            if (i == startNext) {
              return;
            }
            _searchSpace.emplace(containerOption, cellSizeFactor, traversalOption, dataLayoutOption, newton3Option);
            ++i;
          }
        }
      }
    }
  }
}

bool FullSearchMPI::tune(bool currentInvalid) {
  // Repeat as long as traversals are not applicable or we run out of configs
  ++_tuningConfig;

  if (currentInvalid) {
    if (_tuningConfig == _searchSpace.end() || _allLocalConfigsTested) {
      // Will be called several times in case the optimal configuration becomes invalid
      selectLocalOptimalConfiguration();
    }
    return true;
  }

  if (_tuningConfig == _searchSpace.end()) {
    _allLocalConfigsTested = true;

    if (_optimalConfig.cellSizeFactor == -1) {
      // The optimal configuration has not been set yet
      selectLocalOptimalConfiguration();
    }
  }

  // Wait for the Iallreduce from the last tuning step to finish
  // Make all ranks ready for global tuning simultaneously
  AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
  if (_allGlobalConfigsTested) {
    _optimalConfig =
        _configurationCommunicator.optimizeConfiguration(_autopasMPICommunicator, _optimalConfig, _localOptimalTime);
    AutoPasLog(debug, "Global optimal configuration: {}", _optimalConfig.toString());

    return false;
  }

  AutoPas_MPI_Iallreduce(&_allLocalConfigsTested, &_allGlobalConfigsTested, 1, AUTOPAS_MPI_CXX_BOOL, AUTOPAS_MPI_LAND,
                         _autopasMPICommunicator, &_request);

  return true;
}

void FullSearchMPI::selectLocalOptimalConfiguration() {
  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception(
        "FullSearchMPI: Trying to determine fastest configuration without any measurements! "
        "Either selectLocalOptimalConfiguration was called too early or no applicable configurations were found");
  }

  auto optimum = std::min_element(_traversalTimes.begin(), _traversalTimes.end(),
                                  [](std::pair<Configuration, size_t> a, std::pair<Configuration, size_t> b) -> bool {
                                    return a.second < b.second;
                                  });

  AutoPasLog(debug, "Local optimal configuration: {}; took {} nanoseconds", optimum->first.toString(), optimum->second);
  _optimalConfig = optimum->first;
  _localOptimalTime = optimum->second;

  // If the local optimum is searched again with the same traversal times, it should not return the same configuration,
  // because it became invalid
  _traversalTimes.erase(optimum->first);
}

void FullSearchMPI::removeN3Option(Newton3Option badNewton3Option) {
  for (auto ssIter = _searchSpace.begin(); ssIter != _searchSpace.end();) {
    if (ssIter->newton3 == badNewton3Option) {
      // change current config to the next non-deleted
      // if _tuningConfig is in the search space it is assumed that the process is currently tuning
      if (ssIter == _tuningConfig) {
        ssIter = _searchSpace.erase(ssIter);
        _tuningConfig = ssIter;
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

}  // namespace autopas
