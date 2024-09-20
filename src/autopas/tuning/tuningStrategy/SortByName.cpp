/**
 * @file SortByName.cpp
 * @author F. Gratl
 * @date 04.07.23
 */

#include "SortByName.h"
bool autopas::SortByName::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                              const EvidenceCollection &evidenceCollection) {
  // sort configurations to minimize container conversion overhead.
  std::sort(configQueue.begin(), configQueue.end());

  // SortByName does no intentional config wipes to stop the tuning phase
  return false;
}
bool autopas::SortByName::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                const autopas::EvidenceCollection &evidenceCollection) {
  return optimizeSuggestions(configQueue, evidenceCollection);
}

autopas::TuningStrategyOption autopas::SortByName::getOptionType() const { return TuningStrategyOption::sortByName; }
