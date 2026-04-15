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
  // Bump counters and check if they were tuning
  for (const auto &tuner : _autoTuners | std::views::values) {
    tuner->bumpIterationCounters();
  }
}

bool TunerManager::needsLiveInfo(size_t currentIteration) {
  return std::ranges::any_of(_autoTuners, [](const auto &tuner) { return tuner.second->needsLiveInfo(); }) and
         currentIteration < _lastTuningIteration;
}

bool TunerManager::tune(const size_t currentIteration, const LiveInfo &info) {
  if (_lastTuningIteration != currentIteration) {
    _lastTuningIteration = currentIteration;

    for (const auto &tuner : _autoTuners | std::ranges::views::values) {
      if (tuner->needsLiveInfo()) {
        tuner->receiveLiveInfo(info);
      }
    }
    tuneConfigurations();
  }
  return not _tuningFinished;
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

bool TunerManager::tuningPhaseJustFinished() const { return _tuningFinished and (not _tuningJustFinished); }

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

void TunerManager::applyContainerConstraint(std::optional<ContainerOption> containerOption) {
  for (const auto &tuner : _autoTuners | std::views::values) {
    tuner->setContainerConstraint(containerOption);
    // tuner->forceRetune();
    tuner->tuneConfiguration();
  }
}

void TunerManager::tuneConfigurations() {
  // Check if any tuner is still tuning
  bool allTunersFinished = true;
  for (const auto &tuner : _autoTuners | std::views::values) {
    const bool stillTuning = tuner->tuneConfiguration();
    allTunersFinished = allTunersFinished and (not stillTuning);
  }
  _tuningJustFinished = false;
  if (allTunersFinished and (not _tuningFinished)) {
    // Save the best container results before changing container
    setOptimalConfigurations();
    _tuningJustFinished = true;
  }
  _tuningFinished = allTunersFinished;
}

void TunerManager::setOptimalConfigurations() {
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
    AutoPasLog(DEBUG,
               "TunerManager::setOptimalConfigurations: Each AutoTuner will run their individually best rebuild + "
               "traversal configuration");
    applyContainerConstraint(std::nullopt);
  } else {
    AutoPasLog(DEBUG,
               "TunerManager::setOptimalConfigurations: AutoTuners will run their best configuration with the {} "
               "container option",
               bestContainer);
    applyContainerConstraint(bestContainer);
  }
}
}  // namespace autopas