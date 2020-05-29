/**
 * @file FullSearchMPI.h
 * @author W. Thieme
 * @date 17.04.2020
 */

#pragma once

#include <set>
#include <sstream>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Exhaustive full search of the search space by testing every applicable configuration and then selecting the optimum.
 * The search space is distributed among the ranks during the instantiation of a FullSearchMPI object.
 * Every call of tune(false) is globally blocking (i.e. all ranks synchronize).
 * The ranks first search through their local search spaces and then compare their best configurations globally to
 * find the best one in the last step of each tuning phase.
 * Use for homogeneous domains.
 */
class FullSearchMPI : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the FullSearch that generates the rank's portion of the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param comm: provide default value mainly for testing, otherwise TuningStrategyFactory handles this
   */
  FullSearchMPI(const std::set<ContainerOption> &allowedContainerOptions,
                const std::set<double> &allowedCellSizeFactors,
                const std::set<TraversalOption> &allowedTraversalOptions,
                const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                const std::set<Newton3Option> &allowedNewton3Options, AutoPas_MPI_Comm comm = AUTOPAS_MPI_COMM_WORLD);

  /**
   * Constructor for the FullSearch that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   * @param comm: provide default value mainly for testing, otherwise TuningStrategyFactory handles this
   */
  explicit FullSearchMPI(std::set<Configuration> allowedConfigurations, AutoPas_MPI_Comm comm = AUTOPAS_MPI_COMM_WORLD)
      : _containerOptions{},
        _searchSpace(std::move(allowedConfigurations)),
        _tuningConfig(_searchSpace.begin()),
        _optimalConfig(Configuration()),
        _configurationCommunicator(AutoPasConfigurationCommunicator()),
        _localOptimalTime(0),
        _request(AUTOPAS_MPI_REQUEST_NULL),
        _allLocalConfigsTested(false),
        _allGlobalConfigsTested(false),
        _autopasMPICommunicator(comm) {
    for (auto config : _searchSpace) {
      _containerOptions.insert(config.container);
    }
  }

  ~FullSearchMPI() override {
    int finalized;
    AutoPas_MPI_Finalized(&finalized);
    // free the request in case it was still active
    if (not finalized) {
      AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
    }
  }

  inline const Configuration &getCurrentConfiguration() const override {
    // the cellSizeFactor is -1 if _optimalConfig has been initialized to an invalid value
    if (_optimalConfig.cellSizeFactor == -1.) {
      return *_tuningConfig;
    } else {
      return _optimalConfig;
    }
  }

  inline void removeN3Option(Newton3Option badNewton3Option) override;

  inline void addEvidence(long time, size_t /*iteration*/) override { _traversalTimes[*_tuningConfig] = time; }

  inline void reset(size_t /*iteration*/) override {
    _traversalTimes.clear();
    _tuningConfig = _searchSpace.begin();
    _optimalConfig = Configuration();
    _localOptimalTime = 0;
    // set _request to MPI_REQUEST_NULL if it wasn't already
    AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
    _allLocalConfigsTested = false;
    _allGlobalConfigsTested = false;
  }

  inline bool tune(bool = false) override;

  inline std::set<ContainerOption> getAllowedContainerOptions() const override { return _containerOptions; }

  inline bool searchSpaceIsTrivial() const override { return _searchSpace.size() == 1; }

  inline bool searchSpaceIsEmpty() const override { return _searchSpace.empty(); }

 private:
  /**
   * Fills the search space with a subset of the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  inline void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                                  const std::set<double> &allowedCellSizeFactors,
                                  const std::set<TraversalOption> &allowedTraversalOptions,
                                  const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                  const std::set<Newton3Option> &allowedNewton3Options, const int start,
                                  const int startNext);

  /**
   * Uses traversal times from this rank to select the locally optimal configuration.
   */
  inline void selectLocalOptimalConfiguration();

  std::set<ContainerOption> _containerOptions;
  std::set<Configuration> _searchSpace;
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
  std::set<Configuration>::iterator _tuningConfig;
  Configuration _optimalConfig;
  AutoPasConfigurationCommunicator _configurationCommunicator;
  size_t _localOptimalTime;
  AutoPas_MPI_Request _request;
  bool _allLocalConfigsTested;
  bool _allGlobalConfigsTested;
  AutoPas_MPI_Comm _autopasMPICommunicator;
};
}  // namespace autopas
