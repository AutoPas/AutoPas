/**
 * @file SlowConfigFilter.cpp
 * @author F. Gratl
 * @date 29.06.23
 */

#include "SlowConfigFilter.h"

#include <algorithm>

void autopas::SlowConfigFilter::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                    const EvidenceCollection &evidenceCollection) {
  const auto [bestConfig, bestEvidence] = evidenceCollection.getLatestOptimalConfiguration();
  const auto blacklistThreshold = static_cast<long>(static_cast<double>(bestEvidence.value) * _relativeBlacklistRange);
  // push stuff that is slow in the latest tuning phase to the blacklist
  for (const auto &conf : configQueue) {
    const auto *evidenceVec = evidenceCollection.getEvidence(conf);
    // make sure there actually is evidence
    if (evidenceVec and not evidenceVec->empty()) {
      const auto &latestEvidence = evidenceVec->back();
      if (latestEvidence.tuningPhase == bestEvidence.tuningPhase and latestEvidence.value > blacklistThreshold) {
        _blacklist.emplace(conf);
      }
    }
  }

  // Remove all configurations that are on the blacklist
  configQueue.erase(std::remove_if(configQueue.begin(), configQueue.end(),
                                   [&](const auto &conf) { return _blacklist.count(conf) > 0; }),
                    configQueue.end());
}
void autopas::SlowConfigFilter::reset(size_t /*iteration*/, size_t /*tuningPhase*/,
                                      std::vector<Configuration> &configQueue,
                                      const autopas::EvidenceCollection &evidenceCollection) {
  optimizeSuggestions(configQueue, evidenceCollection);
}

autopas::TuningStrategyOption autopas::SlowConfigFilter::getOptionType() {
  return TuningStrategyOption::slowConfigFilter;
}
