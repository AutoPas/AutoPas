/**
 * @file MPIParallelizedStrategy.cpp
 * @author W. Thieme
 * @date 17.09.2020
 */

#include "MPIParallelizedStrategy.h"

namespace autopas {

bool MPIParallelizedStrategy::tune(bool currentInvalid) {
  if (_bucket == AUTOPAS_MPI_COMM_NULL) {
    AutoPasLog(WARN, "_bucket was AUTOPAS_MPI_COMM_NULL");
    _bucket = _comm;
  }

  if (not _strategyStillWorking and currentInvalid) {
    nextFallbackConfig();
    return true;
  }

  if (not _allLocalConfigurationsTested) {
    try {
      _allLocalConfigurationsTested = not _tuningStrategy->tune(currentInvalid);
    } catch (utils::ExceptionHandler::AutoPasException &exception) {
      AutoPasLog(WARN,
                 "MPIParallelizedStrategy: Underlying strategy failed (with error: {}). Reverting to fallback-mode.",
                 exception.what());
      setupFallbackOptions();
    }
  } else if (currentInvalid) {
    AutoPasLog(WARN, "MPIParallelizedStrategy: Underlying strategy found invalid optimum. Reverting to fallback-mode.");
    setupFallbackOptions();
  }

  if (currentInvalid) {
    return true;
  }

  // Wait for the Iallreduce from the last tuning step to finish.
  // Make all ranks ready for global tuning simultaneously.
  // The first time this function is called, it is called with MPI_REQUEST_NULL, where MPI_Wait returns immediately.
  AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
  if (_allGlobalConfigurationsTested) {
    Configuration config = Configuration();
    size_t localOptimalTime = std::numeric_limits<size_t>::max();
    if (_strategyStillWorking) {
      config = _tuningStrategy->getCurrentConfiguration();
      localOptimalTime = _tuningStrategy->getEvidence(config);
    }
    _optimalConfiguration =
        utils::AutoPasConfigurationCommunicator::optimizeConfiguration(_bucket, config, localOptimalTime);

    return false;
  }

  AutoPas_MPI_Iallreduce(&_allLocalConfigurationsTested, &_allGlobalConfigurationsTested, 1, AUTOPAS_MPI_CXX_BOOL,
                         AUTOPAS_MPI_LAND, _bucket, &_request);

  return true;
}

void MPIParallelizedStrategy::setupFallbackOptions() {
  // There was probably an issue finding the optimal configuration.
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
  // Use commSize 1, because then each call advances the configuration exactly by one.
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