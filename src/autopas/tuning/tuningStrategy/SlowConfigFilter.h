/**
 * @file SlowConfigFilter.h
 * @author F. Gratl
 * @date 29.06.23
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "tuning/Configuration.h"
namespace autopas {
/**
 * Acts as a blacklist for configurations that have proven to be very slow.
 */
class SlowConfigFilter : TuningStrategyInterface {
 public:
  explicit SlowConfigFilter(double relativeBlacklistRange) : _relativeBlacklistRange(relativeBlacklistRange){};

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;
  void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;

 private:
  std::set<Configuration> _blacklist{};
  /**
   * Any configuration that is slower than the fastest times this factor will be blacklisted.
   */
  double _relativeBlacklistRange{3};
};

}  // namespace autopas