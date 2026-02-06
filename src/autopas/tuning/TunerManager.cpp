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

void TunerManager::bumpTunerCounters() {
  ++_iteration;
  const bool newTuningPhaseStart = _iteration % _tuningInterval == 0;

  for (const auto &autoTuner : _autoTuners | std::views::values) {
    autoTuner->bumpIterationCounters();
    if (newTuningPhaseStart) {
      autoTuner->incrementTuningPhase();
      autoTuner->setTuningState(true);
    }
  }
}

bool TunerManager::requiresRebuilding() {
  return std::ranges::any_of(_autoTuners,
                             [](const auto &tunerEntry) { return tunerEntry.second->willRebuildNeighborLists(); });
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

}  // namespace autopas