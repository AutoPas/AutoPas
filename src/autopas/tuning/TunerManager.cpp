/**
 * @file TunerManager.cpp
 * @author muehlhaeusser
 * @date 06.02.2026
 */

#include "autopas/tuning/TunerManager.h"

#include <ranges>
#include <set>

namespace autopas {

TunerManager::TunerManager(const AutoTunerInfo &autoTunerInfo) : _tuningInterval(autoTunerInfo.tuningInterval) {}

void TunerManager::addAutoTuner(std::unique_ptr<AutoTuner> tuner, const InteractionTypeOption::Value interactionType) {
  _autoTuners[interactionType] = std::move(tuner);
  setCommonContainerOption();
}

void TunerManager::updateAutoTuners() {
  bumpTunerCounters();
  tuneConfigurations();
}

Configuration TunerManager::rejectConfig(const Configuration &configuration, bool indefinitely) {
  try {
    auto [newConfiguration, stillTuning] =
        _autoTuners.at(configuration.interactionType)->rejectConfig(configuration, indefinitely);
    return newConfiguration;

  } catch (utils::ExceptionHandler::AutoPasException &exception) {
    // Rejected config was the only one for this container
    bool couldChangeContainer = rejectCurrentContainer();
    if (not couldChangeContainer) {
      throw;
    }
    return _autoTuners.at(configuration.interactionType)->getCurrentConfig();
  }
}

bool TunerManager::requiresRebuilding() {
  return std::ranges::any_of(_autoTuners,
                             [](const auto &tunerEntry) { return tunerEntry.second->willRebuildNeighborLists(); });
}

bool TunerManager::allSearchSpacesAreTrivial() const {
  for (const auto &autoTuner : _autoTuners | std::views::values) {
    if (not autoTuner->searchSpaceIsTrivial()) {
      return false;
    }
  }
  return true;
}

bool TunerManager::tuningPhaseJustFinished() const { return _tuningJustFinished; }

void TunerManager::forceRetune() {
  for (const auto &autoTuner : _autoTuners | std::views::values) {
    autoTuner->forceRetune();
  }
  _currentContainerIndex = 0;
  applyContainerConstraint(_commonContainerOptions[_currentContainerIndex]);
}

void TunerManager::setCommonContainerOption() {
  if (_autoTuners.empty()) return;

  // Calculate Intersection of Supported Containers for all interaction types
  bool first = true;
  std::set<ContainerOption> intersectionSet;

  for (const auto &tuner : _autoTuners | std::views::values) {
    std::set<ContainerOption> tunerContainers = tuner->getSearchSpaceContainers();

    if (first) {
      intersectionSet = tunerContainers;
      first = false;
    } else {
      std::set<ContainerOption> result;
      std::ranges::set_intersection(intersectionSet, tunerContainers, std::inserter(result, result.begin()));
      intersectionSet = result;
    }
  }

  // Convert to a vector for deterministic ordering
  _commonContainerOptions.assign(intersectionSet.begin(), intersectionSet.end());

  // Initial Setup: Restrict all tuners to the first container immediately
  if (!_commonContainerOptions.empty()) {
    applyContainerConstraint(_commonContainerOptions[0]);
  }
}

void TunerManager::applyContainerConstraint(ContainerOption containerOption) {
  for (const auto &tuner : _autoTuners | std::views::values) {
    tuner->setContainerConstraint(containerOption);
    tuner->forceRetune();
    tuner->tuneConfiguration();
  }
}

void TunerManager::bumpTunerCounters() {
  ++_iteration;
  bool wasStillTuning = false;

  for (const auto &autoTuner : _autoTuners | std::views::values) {
    wasStillTuning = wasStillTuning or autoTuner->inTuningPhase();
    autoTuner->bumpIterationCounters();
  }
  // New tuning phase is starting
  if (_iteration % _tuningInterval == 0 and not wasStillTuning) {
    _currentContainerIndex = 0;
    applyContainerConstraint(_commonContainerOptions[_currentContainerIndex]);
  }
}

void TunerManager::tuneConfigurations() {
  if (_tuningJustFinished) {
    _tuningJustFinished = false;
  }

  // Check if any tuner is still tuning
  bool wasAlreadyOutsideTuningPhase = true;
  bool allTunersJustFinished = true;
  for (const auto &tuner : _autoTuners | std::views::values) {
    const bool wasTuning = tuner->inTuningPhase();
    const bool stillTuning = tuner->tuneConfiguration();
    wasAlreadyOutsideTuningPhase = wasAlreadyOutsideTuningPhase and (not wasTuning);
    allTunersJustFinished = allTunersJustFinished and (not stillTuning);
  }

  if (allTunersJustFinished and not wasAlreadyOutsideTuningPhase) {
    // Save the best container results before changing container
    captureCurrentContainerPerformance();

    // Try advancing to the next container
    if (_currentContainerIndex + 1 < _commonContainerOptions.size()) {
      ++_currentContainerIndex;
      const ContainerOption nextContainerOption = _commonContainerOptions[_currentContainerIndex];
      applyContainerConstraint(nextContainerOption);
      return;
    }

    selectBestContainer();
  }
}

void TunerManager::captureCurrentContainerPerformance() {
  long totalRuntime = 0;
  const ContainerOption currentContainer = _commonContainerOptions[_currentContainerIndex];

  for (const auto &tuner : _autoTuners | std::views::values) {
    auto [optConf, evidence] = tuner->getEvidenceCollection().getLatestOptimalConfiguration();
    totalRuntime = totalRuntime + evidence.value;
  }

  _containerResults[currentContainer] = totalRuntime;
  AutoPasLog(DEBUG, "TunerManager: Container {} finished with total runtime of {} ns", currentContainer.to_string(),
             totalRuntime);
}

void TunerManager::selectBestContainer() {
  // We have tested all containers. Select the best one. The tuners will find their best config for that container.
  ContainerOption bestContainer = _commonContainerOptions[0];
  long minTime = std::numeric_limits<long>::max();
  for (const auto &[container, time] : _containerResults) {
    if (time < minTime) {
      minTime = time;
      bestContainer = container;
    }
  }

  AutoPasLog(DEBUG, "TunerManager: Tuning Finished. Best Container is {} with total time {} ns.",
             bestContainer.to_string(), minTime);

  for (const auto &tuner : _autoTuners | std::views::values) {
    tuner->setContainerConstraint(bestContainer);
    tuner->selectBestConfiguration();
  }
  _tuningJustFinished = true;
}

bool TunerManager::rejectCurrentContainer() {
  // Remove container from tuning phases
  _commonContainerOptions.erase(_commonContainerOptions.begin() + _currentContainerIndex);

  if (_commonContainerOptions.empty()) {
    // The only allowed container option was rejected.
    return false;
  }

  if (_currentContainerIndex < _commonContainerOptions.size()) {
    // Advance to the next container
    applyContainerConstraint(_commonContainerOptions[_currentContainerIndex]);
  } else {
    // Last container was rejected, so tuning is finished
    selectBestContainer();
  }
  return true;
}

}  // namespace autopas