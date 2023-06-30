/**
 * @file RandomSearch.h
 * @author Jan Nguyen
 * @date 10.07.2019
 */

#pragma once

#include <random>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

/**
 * Randomly test a given number of configurations and select the fastest.
 */
class RandomSearch final : public TuningStrategyInterface {
 public:
  /**
   * Constructor
   * @param maxEvidence Maximum number of evidence this strategy will collect.
   * This number should be smaller than the available search space.
   * @param seed Seed for the random number generator (should only be fixed for deterministic tests).
   */
  explicit RandomSearch(size_t maxEvidence, unsigned long seed);

  void reset(size_t, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Random engine for selecting configurations.
   */
  std::mt19937 _rng;
  /**
   * Since this strategy does not converge this variable defines an upper limit of evidence to collect.
   */
  size_t _maxEvidence;
  /**
   * Counter how many evidence were already collected in this tuning phase.
   */
  size_t _numEvidenceCollected;
};

}  // namespace autopas
