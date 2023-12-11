/**
 * @file SortByName.cpp
 * @author F. Gratl
 * @date 04.07.23
 */

#include "SortByName.h"

void autopas::SortByName::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                              const autopas::EvidenceCollection &evidenceCollection) {
  // sort configurations to minimize container conversion overhead.
  std::sort(configQueue.begin(), configQueue.end());
}
void autopas::SortByName::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                const autopas::EvidenceCollection &evidenceCollection) {
  optimizeSuggestions(configQueue, evidenceCollection);
}

autopas::TuningStrategyOption autopas::SortByName::getOptionType() { return TuningStrategyOption::sortByName; }
