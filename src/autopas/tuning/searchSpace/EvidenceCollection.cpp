/**
 * @file EvidenceCollection.cpp
 * @author F. Gratl
 * @date 28.06.23
 */

#include "EvidenceCollection.h"

#include <limits>

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

std::tuple<Configuration, Evidence> EvidenceCollection::getOptimalConfiguration(std::size_t tuningPhase) const {
  if (_evidenceMap.empty()) {
    utils::ExceptionHandler::exception(
        "EvidenceCollection::getLatestOptimalConfiguration(): Trying to determine the optimal configuration but there "
        "is no evidence yet!");
  }
  Configuration optimalConf{};
  Evidence optimalEvidence{0, 0, std::numeric_limits<decltype(Evidence::value)>::max()};
  for (const auto &[conf, evidenceVec] : _evidenceMap) {
    // reverse iteration of the evidence vector because we are probably interested in the latest evidence.
    for (auto evidenceIter = evidenceVec.rbegin(); evidenceIter != evidenceVec.rend(); ++evidenceIter) {
      if (evidenceIter->tuningPhase == tuningPhase and optimalEvidence.value > evidenceIter->value) {
        optimalConf = conf;
        optimalEvidence = *evidenceIter;
        // Assumption: There is only one evidence per tuning phase.
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