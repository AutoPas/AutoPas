/**
 * @file RandomSearch.h
 * @author F. Gratl
 * @date 17.11.2022
 */

#include "RandomSearch.h"

#include <algorithm>
#include <optional>

#include "autopas/tuning/Configuration.h"

autopas::RandomSearch::RandomSearch(size_t maxEvidence, unsigned long seed)
    : _rng(seed), _maxEvidence(maxEvidence), _numEvidenceCollected(0) {}

void autopas::RandomSearch::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                const EvidenceCollection &evidenceCollection) {
  // Check if enough evidence was collected.
  if (_numEvidenceCollected == _maxEvidence) {
    // Remove everything from the queue and only insert the optimum.
    configQueue.clear();
  } else {
    // Swap a randomly selected configuration to the end.
    std::uniform_int_distribution<std::mt19937::result_type> distribution(0, configQueue.size() - 1);
    const auto selectedIndex = distribution(_rng);
    std::swap(configQueue.back(), configQueue[selectedIndex]);
  }
}

void autopas::RandomSearch::reset(size_t, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                  const autopas::EvidenceCollection &evidenceCollection) {
  _numEvidenceCollected = 0;
}
