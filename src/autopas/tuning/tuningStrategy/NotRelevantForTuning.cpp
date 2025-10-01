/**
 * @file NotRelevantForTuning.cpp
 * @author muehlhaeusser
 * @date 01.10.2025
 */

#include "NotRelevantForTuning.h"

#include "autopas/tuning/Configuration.h"

autopas::NotRelevantForTuning::NotRelevantForTuning() : _numEvidenceCollected(0) {}

bool autopas::NotRelevantForTuning::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                        const EvidenceCollection &evidenceCollection) {
  // Remove everything when a valid config has been found.
  if (_numEvidenceCollected >= 1) {
    configQueue.clear();
    return true;
  }
  // Full Search until a valid config has been found.
  return false;
}

bool autopas::NotRelevantForTuning::reset(size_t, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                          const EvidenceCollection &evidenceCollection) {
  return optimizeSuggestions(configQueue, evidenceCollection);
}

void autopas::NotRelevantForTuning::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  ++_numEvidenceCollected;
};

autopas::TuningStrategyOption autopas::NotRelevantForTuning::getOptionType() const {
  return TuningStrategyOption::notRelevantForTuning;
}
