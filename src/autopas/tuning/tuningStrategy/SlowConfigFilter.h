/**
 * @file SlowConfigFilter.h
 * @author F. Gratl
 * @date 29.06.23
 */

#pragma once

#include <set>
#include <vector>

#include "autopas/options/TuningStrategyOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {
/**
 * Acts as a blacklist for configurations that have proven to be very slow.
 */
class SlowConfigFilter : public TuningStrategyInterface {
 public:
  /**
   * Constructor
   * @param relativeBlacklistRange
   */
  explicit SlowConfigFilter(double relativeBlacklistRange) : _relativeBlacklistRange(relativeBlacklistRange){};

  TuningStrategyOption getOptionType() override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;
  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Set of configurations that should never be used again.
   */
  std::set<Configuration> _blacklist{};
  /**
   * Any configuration that is slower than the fastest times this factor will be blacklisted.
   */
  double _relativeBlacklistRange;
};

}  // namespace autopas