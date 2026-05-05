/**
 * @file TuningManager.cpp
 * @author muehlhaeusser
 * @date 06.02.2026
 */

#include "autopas/tuning/TuningManager.h"

#include <algorithm>
#include <iterator>
#include <ranges>

namespace autopas {

TuningManager::TuningManager(const AutoTunerInfo &autoTunerInfo) : _tuningInterval(autoTunerInfo.tuningInterval) {}

void TuningManager::addAutoTuner(std::unique_ptr<AutoTuner> tuner, const InteractionTypeOption::Value interactionType) {
  _autoTuners[interactionType] = std::move(tuner);
}

bool TuningManager::tune(const size_t currentIteration, const LiveInfo &info) {
  if (_lastTuningIteration != currentIteration) {
    _lastTuningIteration = currentIteration;

    const bool isStart = isStartOfTuningPhase(currentIteration);
    const bool aboutToBegin = tuningPhaseAboutToBegin(currentIteration);
    if (isStart) {
      ++_tuningPhase;
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
  if (_autoTuners.at(interactionType)->inTuningPhase()) {
    _autoTuners.at(interactionType)
        ->addMeasurement(sampleRebuild, sampleTraverseParticles, neighborListRebuilt, iteration, _tuningPhase);
  }
}

void TuningManager::setRebuildFrequency(const double rebuildFrequency) {
  for (const auto &tuner : _autoTuners | std::views::values) {
    tuner->setRebuildFrequency(rebuildFrequency);
  }
}

Configuration TuningManager::rejectConfiguration(const Configuration &rejectedConfig, const bool indefinitely,
                                                 const InteractionTypeOption::Value interactionType) {
  auto &tunerToReject = _autoTuners[interactionType];
  auto newConfig = tunerToReject->rejectConfig(rejectedConfig, indefinitely, _tuningPhase);
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
  const bool needLiveInfo =
      isStart and std::ranges::any_of(_autoTuners, [&](const auto &tuner) { return tuner.second->needsLiveInfo(); });
  const bool needDomainStats = aboutToBegin and std::ranges::any_of(_autoTuners, [&](const auto &tuner) {
                                 return tuner.second->needsDomainSimilarityStatistics();
                               });
  return needLiveInfo or needDomainStats;
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

bool TuningManager::isStartOfTuningPhase(const size_t currentIteration) const {
  // If it's a multiple of the interval, or if a retune was forced
  return (currentIteration % _tuningInterval == 0) or _forceRetunePending;
}

bool TuningManager::inFirstConfigurationLastSample(const size_t currentIteration) const {
  return (currentIteration % _tuningInterval == _maxSamples - 1);
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
  // 1. Check the potential time if we switch container with every functor call
  long bestTotalValue = 0;
  for (const auto &tuner : _autoTuners | std::views::values) {
    auto [optConf, evidence] = tuner->getEvidenceCollection().getBestConfigWithRebuild();
    bestTotalValue += evidence.rebuildValue + evidence.traversalValue;
  }

  // 2. Check the potential time if we keep all interaction types fixed to the same container (and save rebuilds)
  long bestCommonContainerValue = std::numeric_limits<long>::max();
  auto commonContainerAndCellSizeFactorOptions = getCommonContainerAndCellSizeFactors();

  std::optional<ContainerOption> bestContainer = std::nullopt;
  std::optional<double> bestCSF = std::nullopt;

  for (const auto &[container, cellSizeFactor] : commonContainerAndCellSizeFactorOptions) {
    long totalValueForContainerAndCSF = 0;
    bool containerWithCSFIsValid = true;

    for (const auto &tuner : _autoTuners | std::views::values) {
      try {
        auto [optConf, evidence] =
            tuner->getEvidenceCollection().getBestConfigForContainerAndCSF(container, cellSizeFactor);
        totalValueForContainerAndCSF += evidence.effectiveValue;
      } catch (const std::exception &e) {
        containerWithCSFIsValid = false;
        AutoPasLog(DEBUG, "No evidence for the interaction type {} with container {} and cell size factor {}",
                   interactionType, container, cellSizeFactor);
        break;
      }
    }
    if (containerWithCSFIsValid and totalValueForContainerAndCSF < bestCommonContainerValue) {
      bestCommonContainerValue = totalValueForContainerAndCSF;
      bestContainer = container;
      bestCSF = cellSizeFactor;
    }
  }

  // 3. Evaluate our options
  if ((not bestContainer.has_value()) or (not bestCSF.has_value()) or bestTotalValue < bestCommonContainerValue) {
    AutoPasLog(DEBUG,
               "TuningManager::setOptimalConfigurations: Each AutoTuner will run their individually best rebuild + "
               "traversal configuration");
    for (const auto &tuner : _autoTuners | std::views::values) {
      auto [optConf, evidence] = tuner->getEvidenceCollection().getBestConfigWithRebuild();
      tuner->setOptimalConfiguration(optConf);
    }
  } else {
    AutoPasLog(DEBUG,
               "TuningManager::setOptimalConfigurations: AutoTuners will run their best configuration with the {} "
               "container option",
               bestContainer.value());
    for (const auto &tuner : _autoTuners | std::views::values) {
      auto [optConf, evidence] =
          tuner->getEvidenceCollection().getBestConfigForContainerAndCSF(bestContainer.value(), bestCSF.value());
      tuner->setOptimalConfiguration(optConf);
    }
  }
}

std::set<std::pair<ContainerOption, double>> TuningManager::getCommonContainerAndCellSizeFactors() const {
  if (_autoTuners.empty()) return {};

  // Calculate Intersection of Supported Containers for all interaction types
  std::set<ContainerOption> commonContainers;
  std::set<double> commonCellSizeFactors;
  bool first = true;
  std::set<std::pair<ContainerOption, double>> intersectionSet;

  for (const auto &tuner : _autoTuners | std::views::values) {
    std::set<ContainerOption> tunerContainers = tuner->getSearchSpaceContainers();
    std::set<double> tunerCellSizeFactors = tuner->getSearchSpaceCellSizeFactors();

    if (first) {
      commonContainers = std::move(tunerContainers);
      commonCellSizeFactors = std::move(tunerCellSizeFactors);
      first = false;
    } else {
      // Intersect Container Options
      std::set<ContainerOption> resultContainers;
      std::ranges::set_intersection(commonContainers, tunerContainers,
                                    std::inserter(resultContainers, resultContainers.begin()));
      commonContainers = std::move(resultContainers);

      // Intersect Cell Size Factors
      std::set<double> resultFactors;
      std::ranges::set_intersection(commonCellSizeFactors, tunerCellSizeFactors,
                                    std::inserter(resultFactors, resultFactors.begin()));
      commonCellSizeFactors = std::move(resultFactors);
    }

    // Early exit if either common set becomes empty
    if (commonContainers.empty() || commonCellSizeFactors.empty()) {
      return {};
    }
  }

  // Generate the combinations (Cartesian product) of the remaining valid elements
  std::set<std::pair<ContainerOption, double>> finalCombinations;
  for (const auto &container : commonContainers) {
    for (const auto &factor : commonCellSizeFactors) {
      finalCombinations.emplace(container, factor);
    }
  }

  return finalCombinations;
}

bool TuningManager::tuningPhaseAboutToBegin(const size_t currentIteration) const {
  return _tuningInterval <= 10 or currentIteration % _tuningInterval >= _tuningInterval - 10;
}
}  // namespace autopas