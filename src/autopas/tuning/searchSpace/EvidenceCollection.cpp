/**
 * @file EvidenceCollection.cpp
 * @author F. Gratl
 * @date 28.06.23
 */

#include "EvidenceCollection.h"

#include <limits>
#include <ranges>

#include "Evidence.h"
#include "autopas/tuning/Configuration.h"

namespace autopas {

void EvidenceCollection::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  _latestTuningPhase = std::max(_latestTuningPhase, evidence.tuningPhase);
  _evidenceMap[configuration].push_back(evidence);
}

const std::vector<Evidence> *EvidenceCollection::getEvidence(const Configuration &configuration) const {
  const auto iter = _evidenceMap.find(configuration);
  if (iter != _evidenceMap.end()) {
    return &iter->second;
  }
  return nullptr;
}

Evidence &EvidenceCollection::modifyLastEvidence(const Configuration &configuration) {
  return _evidenceMap[configuration].back();
}

std::tuple<Configuration, Evidence> EvidenceCollection::getBestConfigForContainer(
    ContainerOption containerOption) const {
  return getOptimalConfiguration(_latestTuningPhase, EvidenceMode::REDUCED, containerOption);
}

std::tuple<Configuration, Evidence> EvidenceCollection::getBestConfigNotReduced() const {
  return getOptimalConfiguration(_latestTuningPhase, EvidenceMode::TOTAL, std::nullopt);
}

std::tuple<Configuration, Evidence> EvidenceCollection::getOptimalConfiguration(
    const size_t tuningPhase, const EvidenceMode mode, const std::optional<ContainerOption> containerConstraint) const {
  if (_evidenceMap.empty()) {
    utils::ExceptionHandler::exception(
        "EvidenceCollection::getOptimalConfiguration(): Trying to determine the optimal configuration but there "
        "is no evidence yet!");
  }
  Configuration optimalConf{};
  Evidence optimalEvidence{};
  long bestValue = std::numeric_limits<long>::max();

  // Helper lambda to fetch the metric as determined by the EvidenceMode
  auto getValue = [mode](const Evidence &e) -> long {
    switch (mode) {
      case EvidenceMode::REDUCED:
        return e.reducedValue;
      case EvidenceMode::TRAVERSAL:
        return e.traversalValue;
      case EvidenceMode::TOTAL:
        return e.rebuildValue + e.traversalValue;
    }
    return std::numeric_limits<long>::max();
  };

  for (const auto &[conf, evidenceVec] : _evidenceMap) {
    if (containerConstraint.has_value() and conf.container != containerConstraint.value()) {
      continue;
    }
    // reverse iteration of the evidence vector because we are probably interested in the latest evidence.
    for (const auto &evidence : std::views::reverse(evidenceVec)) {
      if (evidence.tuningPhase == tuningPhase) {
        long currentValue = getValue(evidence);

        if (currentValue < bestValue) {
          bestValue = currentValue;
          optimalConf = conf;
          optimalEvidence = evidence;
        }

        // Assumption: There is only one evidence per tuning phase
        break;
      }
    }
  }
  // Check if the found config differs from the (invalid) default. If not nothing was found.
  if (optimalConf == Configuration{}) {
    utils::ExceptionHandler::exception(
        "EvidenceCollection::getLatestOptimalConfiguration(): No configuration could be determined to be the optimum "
        "for tuning phase {}. This suggests there is no evidence for this phase yet.",
        tuningPhase);
  }
  return {optimalConf, optimalEvidence};
}

std::tuple<Configuration, Evidence> EvidenceCollection::getLatestOptimalConfiguration() const {
  return getOptimalConfiguration(_latestTuningPhase);
}

bool EvidenceCollection::empty() const { return _evidenceMap.empty(); }
}  // namespace autopas