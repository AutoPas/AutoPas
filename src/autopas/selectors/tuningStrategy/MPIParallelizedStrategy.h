/**
 * @file MPIParallelizedStrategy.h
 * @author W. Thieme
 * @date 27.05.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "autopas/utils/ConfigurationAndRankIteratorHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas {

/**
 * Wrapper for other tuning strategies which handles MPI communication between processes.
 * The splitting of the search space is not currently handled by this class, but by AutoPasConfigurationCommunicator.
 * This mainly overwrites the tune() method to globally compare the best configuration.
 */
class MPIParallelizedStrategy : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the wrapper. Assumes that the tuningStrategy has already been constructed with the appropriate
   * search space.
   * If the underlying tuning strategy fails for some reason, uses the fallback-options to keep global search going.
   * @param tuningStrategy The underlying tuning strategy which tunes locally on it's node.
   * @param comm The communicator holding all ranks which participate in this tuning strategy
   * @param fallbackContainers
   * @param fallbackTraversals
   * @param fallbackLoadEstimators
   * @param fallbackDataLayouts
   * @param fallbackNewton3s
   */
  MPIParallelizedStrategy(std::unique_ptr<TuningStrategyInterface> tuningStrategy, const AutoPas_MPI_Comm comm,
                          const std::set<autopas::ContainerOption> &fallbackContainers,
                          const std::set<autopas::TraversalOption> &fallbackTraversals,
                          const std::set<autopas::LoadEstimatorOption> &fallbackLoadEstimators,
                          const std::set<autopas::DataLayoutOption> &fallbackDataLayouts,
                          const std::set<autopas::Newton3Option> &fallbackNewton3s)
      : _tuningStrategy(std::move(tuningStrategy)),
        _comm(comm),
        _fallbackContainers(fallbackContainers),
        _fallbackTraversalOptions(fallbackTraversals),
        _fallbackLoadEstimators(fallbackLoadEstimators),
        _fallbackDataLayouts(fallbackDataLayouts),
        _fallbackNewton3s(fallbackNewton3s) {}

  inline void addEvidence(long time, size_t iteration) override { _tuningStrategy->addEvidence(time, iteration); }

  [[nodiscard]] inline long getEvidence(Configuration configuration) const override {
    return _tuningStrategy->getEvidence(configuration);
  }

  [[nodiscard]] const Configuration &getCurrentConfiguration() const override {
    // cellSizeFactor == -1 iff a config is invalid
    if (_optimalConfiguration.cellSizeFactor == -1) {
      return _tuningStrategy->getCurrentConfiguration();
    } else {
      return _optimalConfiguration;
    }
  }

  bool tune(bool currentInvalid) override;

  void reset(size_t iteration) override {
    _optimalConfiguration = Configuration();
    _allGlobalConfigurationsTested = false;
    _allLocalConfigurationsTested = false;
    _strategyStillWorking = true;
    if (_configIterator != nullptr) {
      _configIterator.reset();
    }
    try {
      _tuningStrategy->reset(iteration);
    } catch (utils::ExceptionHandler::AutoPasException &exception) {
      AutoPasLog(warn, "MPIParallelizedStrategy: Underlying strategy failed. Reverting to fallback-mode.");
      setupFallbackOptions();
    }
  }

  [[nodiscard]] std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _tuningStrategy->getAllowedContainerOptions();
  }

  inline void removeN3Option(Newton3Option badN3Option) override { _tuningStrategy->removeN3Option(badN3Option); }

  [[nodiscard]] inline bool searchSpaceIsTrivial() const override { return _tuningStrategy->searchSpaceIsTrivial(); }

  [[nodiscard]] inline bool searchSpaceIsEmpty() const override { return _tuningStrategy->searchSpaceIsEmpty(); }

  /**
   * Getter for the internal tuningStrategy
   * @return the tuning strategy which was provided in the constructor
   */
  [[nodiscard]] inline const TuningStrategyInterface &getTuningStrategy() const { return *_tuningStrategy; }

 private:
  /**
   * Essentially sets up a full-search-esque configuration selection, in case the underlying search strategy fails.
   * Usually the strategy fails, because it has no valid configurations when many ranks are present.
   */
  void setupFallbackOptions();

  /**
   * Uses _configIterator to get a new config from the global search space.
   */
  void nextFallbackConfig();

  /**
   * The tuning strategy tuning locally
   */
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  AutoPas_MPI_Comm _comm;
  AutoPas_MPI_Request _request{AUTOPAS_MPI_REQUEST_NULL};

  /**
   * The globally optimal configuration.
   * Usually holds a value that is not in a given rank's search space
   */
  Configuration _optimalConfiguration{Configuration()};

  bool _allLocalConfigurationsTested{false};
  bool _allGlobalConfigurationsTested{false};
  bool _strategyStillWorking{true};

  // fallback configurations, in case the underlying search strategy fails.
  const std::set<ContainerOption> _fallbackContainers;
  // fallback options cannot deal with continuous cellSizeFactors
  const autopas::NumberSetFinite<double> _fallbackCellSizeFactor{1};
  const std::set<TraversalOption> _fallbackTraversalOptions;
  const std::set<LoadEstimatorOption> _fallbackLoadEstimators;
  const std::set<DataLayoutOption> _fallbackDataLayouts;
  const std::set<Newton3Option> _fallbackNewton3s;
  int _numFallbackConfigs{-1};
  std::unique_ptr<utils::ConfigurationAndRankIteratorHandler> _configIterator{nullptr};
};

bool MPIParallelizedStrategy::tune(bool currentInvalid) {
  int rank;
  AutoPas_MPI_Comm_rank(_comm, &rank);

  if (not _strategyStillWorking and currentInvalid) {
    nextFallbackConfig();
    return true;
  }

  if (not _allLocalConfigurationsTested) {
    try {
      _allLocalConfigurationsTested = not _tuningStrategy->tune(currentInvalid);
    } catch (utils::ExceptionHandler::AutoPasException &exception) {
      AutoPasLog(warn, "MPIParallelizedStrategy: Underlying strategy failed. Reverting to fallback-mode.");
      setupFallbackOptions();
    }
  } else if (currentInvalid) {
    AutoPasLog(warn, "MPIParallelizedStrategy: Underlying strategy found invalid optimum. Reverting to fallback-mode.");
    setupFallbackOptions();
  }

  if (currentInvalid) {
    return true;
  }

  // Wait for the Iallreduce from the last tuning step to finish
  // Make all ranks ready for global tuning simultaneously
  AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
  if (_allGlobalConfigurationsTested) {
    Configuration config = Configuration();
    size_t localOptimalTime = std::numeric_limits<size_t>::max();
    if (_strategyStillWorking) {
      config = _tuningStrategy->getCurrentConfiguration();
      localOptimalTime = _tuningStrategy->getEvidence(config);
    }
    _optimalConfiguration =
        utils::AutoPasConfigurationCommunicator::optimizeConfiguration(_comm, config, localOptimalTime);

    return false;
  }

  AutoPas_MPI_Iallreduce(&_allLocalConfigurationsTested, &_allGlobalConfigurationsTested, 1, AUTOPAS_MPI_CXX_BOOL,
                         AUTOPAS_MPI_LAND, _comm, &_request);

  return true;
}

void MPIParallelizedStrategy::setupFallbackOptions() {
  // There was probably an issue finding the optimal configuration
  _allLocalConfigurationsTested = true;
  _strategyStillWorking = false;

  // Essentially turn into full search if the underlying strategy dies.
  if (_numFallbackConfigs == -1) {
    _numFallbackConfigs = utils::AutoPasConfigurationCommunicator::getSearchSpaceSize(
        _fallbackContainers, _fallbackCellSizeFactor, _fallbackTraversalOptions, _fallbackLoadEstimators,
        _fallbackDataLayouts, _fallbackNewton3s);
  }
  auto numbersSet = _fallbackCellSizeFactor.getAll();
  _configIterator = std::make_unique<utils::ConfigurationAndRankIteratorHandler>(
      _fallbackContainers, numbersSet, _fallbackTraversalOptions, _fallbackLoadEstimators, _fallbackDataLayouts,
      _fallbackNewton3s, _numFallbackConfigs, 1);
  _optimalConfiguration =
      Configuration(*_configIterator->getContainerIterator(), *_configIterator->getCellSizeFactorIterator(),
                    *_configIterator->getTraversalIterator(), *_configIterator->getLoadEstimatorIterator(),
                    *_configIterator->getDataLayoutIterator(), *_configIterator->getNewton3Iterator());
}

void MPIParallelizedStrategy::nextFallbackConfig() {
  _configIterator->advanceIterators(_numFallbackConfigs, 1);
  if (_configIterator->getRankIterator() >= 1) {
    // All strategies have been searched through and rejected.
    utils::ExceptionHandler::exception("MPIParallelizedStrategy: No viable configurations were provided.");
  }
  _optimalConfiguration =
      Configuration(*_configIterator->getContainerIterator(), *_configIterator->getCellSizeFactorIterator(),
                    *_configIterator->getTraversalIterator(), *_configIterator->getLoadEstimatorIterator(),
                    *_configIterator->getDataLayoutIterator(), *_configIterator->getNewton3Iterator());
}

}  // namespace autopas
