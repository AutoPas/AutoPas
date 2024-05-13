/**
 * @file SortByName.h
 * @author F. Gratl
 * @date 04.07.23
 */

#pragma once

#include "TuningStrategyInterface.h"

namespace autopas {

/**
 * This strategy sorts the given queue by Configuration::operator<() to minimize the container conversion overhead.
 *
 * @note Only use this if you plan to test the full queue, otherwise optimal configurations might be sorted beyond
 * the point you test.
 */
class SortByName : public TuningStrategyInterface {
 public:
  TuningStrategyOption getOptionType() const override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;
  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;
};

}  // namespace autopas
