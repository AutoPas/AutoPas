/**
 * @file TuningManager.cpp
 * @author muehlhaeusser
 * @date 06.02.2026
 */

#include "autopas/tuning/TuningManager.h"

#include <ranges>
#include <set>

namespace autopas {

TuningManager::TuningManager(const AutoTunerInfo &autoTunerInfo) : _tuningInterval(autoTunerInfo.tuningInterval) {}

void TuningManager::addAutoTuner(std::unique_ptr<AutoTuner> tuner, const InteractionTypeOption::Value interactionType) {
  _autoTuners[interactionType] = std::move(tuner);
  setCommonContainerOption();
}

bool TuningManager::tune(const size_t currentIteration, const LiveInfo &info) {
  if (_lastTuningIteration != currentIteration) {
    _lastTuningIteration = currentIteration;

    // Evaluate if we crossed the interval threshold
    if (currentIteration % _tuningInterval == 0) {
      ++_tuningPhase;
    }

    const bool isStart = isStartOfTuningPhase(currentIteration);
    const bool aboutToBegin = tuningPhaseAboutToBegin(currentIteration);
    if (isStart) {
      _forceRetunePending = false;
    }

    for (const auto &tuner : _autoTuners | std::ranges::views::values) {
      if (tuner->needsLiveInfo() and (isStart or aboutToBegin)) {
        tuner->receiveLiveInfo(info, isStart);
      }
    }
    tuneConfigurations(currentIteration);
  }
  return not _tuningFinished;
}

void TuningManager::addMeasurement(long sampleRebuild, long sampleTraverseParticles, bool neighborListRebuilt,
                                   size_t iteration, InteractionTypeOption::Value interactionType) const {
  _autoTuners.at(interactionType)
      ->addMeasurement(sampleRebuild, sampleTraverseParticles, neighborListRebuilt, iteration, _tuningPhase);
}

Configuration TuningManager::rejectConfiguration(const Configuration &rejectedConfig, bool indefinitely,
                                                 size_t currentIteration,
                                                 InteractionTypeOption::Value interactionType) {
  auto &tunerToReject = _autoTuners[interactionType];
  auto newConfig = tunerToReject->rejectConfig(rejectedConfig, indefinitely, currentIteration, _tuningPhase);
  bool allTunersFinished = true;
  for (const auto &tuner : _autoTuners | std::views::values) {
    allTunersFinished = allTunersFinished and (not tuner->inTuningPhase());
  }
  _tuningFinished = allTunersFinished;
  return newConfig;
}

void TuningManager::forceRetune() {
  _forceRetunePending = true;
  for (const auto &autoTuner : _autoTuners | std::views::values) {
    autoTuner->forceRetune();
  }
}

void TuningManager::logTuningResult(long tuningTime, size_t currentIteration,
                                    InteractionTypeOption::Value interactionType) const {
  _autoTuners.at(interactionType)->logTuningResult(tuningTime, currentIteration);
}

bool TuningManager::requiresRebuilding(const size_t currentIteration) {
  const bool isStart = isStartOfTuningPhase(currentIteration);
  return isStart or std::ranges::any_of(_autoTuners, [&](const auto &tunerEntry) {
           return tunerEntry.second->willRebuildNeighborLists();
         });
}

bool TuningManager::needsLiveInfo(const size_t currentIteration) {
  const bool isStart = isStartOfTuningPhase(currentIteration);
  const bool aboutToBegin = tuningPhaseAboutToBegin(currentIteration);
  return (isStart or aboutToBegin) and
         std::ranges::any_of(_autoTuners, [&](const auto &tuner) { return tuner.second->needsLiveInfo(); });
}

bool TuningManager::allSearchSpacesAreTrivial() const {
  for (const auto &autoTuner : _autoTuners | std::views::values) {
    if (not autoTuner->searchSpaceIsTrivial()) {
      return false;
    }
  }
  return true;
}

bool TuningManager::tuningPhaseJustFinished() const {
  return _tuningFinished and (not _transitionToOptimalConfigurations);
}

bool TuningManager::inFirstTuningIteration(const size_t currentIteration) const {
  return currentIteration % _tuningInterval == 0;
}

const Configuration &TuningManager::getCurrentConfig(const InteractionTypeOption::Value interactionType) const {
  return _autoTuners.at(interactionType)->getCurrentConfig();
}

const TuningMetricOption &TuningManager::getTuningMetric(const InteractionTypeOption::Value interactionType) const {
  return _autoTuners.at(interactionType)->getTuningMetric();
}

void TuningManager::tuneConfigurations(const size_t currentIteration) {
  // Check if any tuner is still tuning
  bool allTunersFinished = true;
  const bool isStart = isStartOfTuningPhase(currentIteration);

  for (const auto &tuner : _autoTuners | std::views::values) {
    const bool stillTuning = tuner->tuneConfiguration(currentIteration, _tuningPhase, isStart);
    allTunersFinished = allTunersFinished and (not stillTuning);
  }
  _transitionToOptimalConfigurations = false;
  if (allTunersFinished and (not _tuningFinished)) {
    // Save the best container results before changing container
    setOptimalConfigurations();
    _transitionToOptimalConfigurations = true;
  }
  _tuningFinished = allTunersFinished;
}

void TuningManager::setOptimalConfigurations() {
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
               "TuningManager::setOptimalConfigurations: Each AutoTuner will run their individually best rebuild + "
               "traversal configuration");
    for (const auto &tuner : _autoTuners | std::views::values) {
      auto [optConf, evidence] = tuner->getEvidenceCollection().getBestConfigNotReduced();
      tuner->forceOptimalConfiguration(optConf);
    }
  } else {
    AutoPasLog(DEBUG,
               "TuningManager::setOptimalConfigurations: AutoTuners will run their best configuration with the {} "
               "container option",
               bestContainer);
    for (const auto &tuner : _autoTuners | std::views::values) {
      auto [optConf, evidence] = tuner->getEvidenceCollection().getBestConfigForContainer(bestContainer);
      tuner->forceOptimalConfiguration(optConf);
    }
  }
}

std::set<ContainerOption> TuningManager::setCommonContainerOption() {
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

bool TuningManager::tuningPhaseAboutToBegin(const size_t currentIteration) const {
  return _tuningInterval <= 10 or currentIteration % _tuningInterval > _tuningInterval - 10;
}

bool TuningManager::isStartOfTuningPhase(const size_t currentIteration) const {
  // If it's a multiple of the interval, or if a retune was forced
  return (currentIteration % _tuningInterval == 0) or _forceRetunePending;
}

}  // namespace autopas