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
  } else {
    return nullptr;
  }
}

Evidence &EvidenceCollection::modifyLastEvidence(const Configuration &configuration) {
  return _evidenceMap[configuration].back();
}

std::tuple<Configuration, Evidence> EvidenceCollection::getBestConfigForContainer(
    ContainerOption containerOption) const {
  return getOptimalConfiguration(_latestTuningPhase, EvidenceMode::REDUCED, containerOption);
}

std::tuple<Configuration, Evidence> EvidenceCollection::getBestConfigNotReduced() const {
  if (_evidenceMap.empty()) {
    utils::ExceptionHandler::exception(
        "EvidenceCollection::getBestConfigNotReduced(): Trying to determine the optimal configuration but there "
        "is no evidence yet!");
  }
  Configuration optimalConf{};
  Evidence optimalEvidence{0, 0, std::numeric_limits<decltype(Evidence::reducedValue)>::max(),
                           std::numeric_limits<decltype(Evidence::rebuildValue)>::max() / 2,
                           std::numeric_limits<decltype(Evidence::traversalValue)>::max() / 2};

  for (const auto &[conf, evidenceVec] : _evidenceMap) {
    // reverse iteration of the evidence vector because we are probably interested in the latest evidence.
    for (const auto &evidenceIter : std::views::reverse(evidenceVec)) {
      if (evidenceIter.tuningPhase == _latestTuningPhase and
          (optimalEvidence.rebuildValue + optimalEvidence.traversalValue) >
              (evidenceIter.rebuildValue + evidenceIter.traversalValue)) {
        optimalConf = conf;
        optimalEvidence = evidenceIter;
        // Assumption: There is only one evidence per tuning phase.
        break;
      }
    }
  }
  // Check if the found config differs from the (invalid) default. If not nothing was found.
  if (optimalConf == Configuration{}) {
    utils::ExceptionHandler::exception(
        "EvidenceCollection::getBestConfigNotReduced(): No configuration could be determined to be the optimum "
        "for tuning phase {}. This suggests there is no evidence for this phase yet.",
        _latestTuningPhase);
  }
  return {optimalConf, optimalEvidence};
}

std::tuple<Configuration, Evidence> EvidenceCollection::getOptimalConfiguration(
    size_t tuningPhase, EvidenceMode mode, const std::optional<ContainerOption> containerConstraint) const {
  if (_evidenceMap.empty()) {
    utils::ExceptionHandler::exception(
        "EvidenceCollection::getOptimalConfiguration(): Trying to determine the optimal configuration but there "
        "is no evidence yet!");
  }
  Configuration optimalConf{};
  Evidence optimalEvidence{0, 0, std::numeric_limits<long>::max(), std::numeric_limits<long>::max() / 2,
                           std::numeric_limits<long>::max() / 2};
  for (const auto &[conf, evidenceVec] : _evidenceMap) {
    if (containerConstraint.has_value() and conf.container != containerConstraint.value()) {
      continue;
    }
    // reverse iteration of the evidence vector because we are probably interested in the latest evidence.
    for (const auto &evidence : std::views::reverse(evidenceVec)) {
      if (evidence.tuningPhase == tuningPhase) {
        switch (mode) {
          case EvidenceMode::REDUCED:
            if (optimalEvidence.reducedValue > evidence.reducedValue) {
              optimalConf = conf;
              optimalEvidence = evidence;
            }
            break;
          case EvidenceMode::TRAVERSAL:
            if (optimalEvidence.traversalValue > evidence.traversalValue) {
              optimalConf = conf;
              optimalEvidence = evidence;
            }
            break;
          case EvidenceMode::TOTAL:
            if ((optimalEvidence.rebuildValue + optimalEvidence.traversalValue) >
                (evidence.rebuildValue + evidence.traversalValue)) {
              optimalConf = conf;
              optimalEvidence = evidence;
            }
        }
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