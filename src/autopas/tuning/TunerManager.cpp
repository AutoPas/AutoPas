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

void TunerManager::updateAutoTuners(const size_t currentIteration) {
  bool anyTunerWasTuning = false;

  // Bump counters and check if they were tuning
  for (const auto &tuner : _autoTuners | std::views::values) {
    anyTunerWasTuning = anyTunerWasTuning or tuner->inTuningPhase();
    tuner->bumpIterationCounters();
  }

  // Tune configurations
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
}

std::set<ContainerOption> TunerManager::setCommonContainerOption() {
  if (_autoTuners.empty()) return {};

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

  return intersectionSet;
}

void TunerManager::applyContainerConstraint(ContainerOption containerOption) {
  for (const auto &tuner : _autoTuners | std::views::values) {
    tuner->setContainerConstraint(containerOption);
    tuner->forceRetune();
    tuner->tuneConfiguration();
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
    findBestConfigsCombination();
  }
}

void TunerManager::findBestConfigsCombination() {

  long bestTotalValue = 0;
  for (const auto &tuner : _autoTuners | std::views::values) {
    auto [optConf, evidence] = tuner->getEvidenceCollection().getBestConfigNotReduced();
    bestTotalValue += evidence.rebuildValue + evidence.traversalValue;
  }
  long bestCommonContainerValue = std::numeric_limits<long>::max();
  auto commonContainerOptions = setCommonContainerOption();
  ContainerOption bestContainer = *commonContainerOptions.begin();

  for (const auto &container : commonContainerOptions) {
    long totalValueForContainer = 0;
    for (const auto &tuner : _autoTuners | std::views::values) {
      auto [optConf, evidence] = tuner->getEvidenceCollection().getBestConfigForContainer(container);
      totalValueForContainer += evidence.reducedValue;
    }
    if (totalValueForContainer < bestCommonContainerValue) {
      bestCommonContainerValue = totalValueForContainer;
      bestContainer = container;
    }
  }

  if (bestTotalValue < bestCommonContainerValue) {
    AutoPasLog(DEBUG, "TunerManager::findBestConfigsCombination: ");
    // ContainerConstraint = None
  } else {
    AutoPasLog(DEBUG, "TunerManager::findBestConfigsCombination: ");
    applyContainerConstraint(bestContainer);
  }

}
}  // namespace autopas