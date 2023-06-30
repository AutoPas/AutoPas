/**
 * @file EvidenceCollection.h
 * @author F. Gratl
 * @date 28.06.23
 */

#pragma once

#include <cstddef>
#include <map>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/Evidence.h"

namespace autopas {
/**
 * Class to manage all evidence.
 */
class EvidenceCollection {
 public:
  EvidenceCollection() = default;

  /**
   * Store a piece of evidence in the internal storage.
   * @param configuration
   * @param evidence
   */
  void addEvidence(const Configuration &configuration, const Evidence &evidence);

  /**
   * Returns all evidence collected for a given configuration.
   * @param configuration
   * @return Const reference to the vector of collected evidence.
   */
  const std::vector<Evidence> &getEvidence(const Configuration &configuration) const;

  /**
   * Returns a modifiable reference to the last evidence of a given configuration.
   * @param configuration
   * @return Non-const reference to the last entry in the vector of a given configuration.
   */
  Evidence &modifyLastEvidence(const Configuration &configuration);

  /**
   * Retrieve the configuration with the lowest evidence value for the latest tuning phase.
   *
   * @return The optimal configuration.
   */
  std::tuple<Configuration, Evidence> getLatestOptimalConfiguration() const;

 private:
  /**
   * For each configuration the collection of all evidence (smoothed values) collected so far and in which iteration.
   * Configuration -> vector< iteration, time >
   */
  std::map<Configuration, std::vector<Evidence>> _evidenceMap{};

  /**
   * Number of the newest tuning phase for which evidence exists in this object.
   * This tuning phase is not necessarily completed yet.
   */
  size_t _latestTuningPhase;
};
}  // namespace autopas