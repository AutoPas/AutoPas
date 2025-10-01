/**
 * @file NotRelevantForTuning.h
 * @author muehlhaeusser
 * @date 01.10.2025
 */

#pragma once

#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

/**
 * Look for any valid configuration.
 */
class NotRelevantForTuning final : public TuningStrategyInterface {
 public:
  /**
   * Constructor
   */
  explicit NotRelevantForTuning();

  TuningStrategyOption getOptionType() const override;

  bool reset(size_t, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

 private:
  /**
   * Counter how much evidence was already collected in this tuning phase.
   */
  size_t _numEvidenceCollected;
};

}  // namespace autopas
