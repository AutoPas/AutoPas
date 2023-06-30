/**
 * @file TuningStrategyInterface.h
 * @author F. Gratl
 * @date 29.05.2019
 */

#pragma once

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/LiveInfo.h"

namespace autopas {

/**
 * Interface for tuning strategies for the auto tuner.
 */
class TuningStrategyInterface {
 public:
  virtual ~TuningStrategyInterface() = default;

  /**
   * Store empirically collected information for the current configuration.
   *
   * Implementing this function is only necessary if the tuning strategy processes evidence differently
   * than EvidenceCollection.
   *
   * @param configuration Measured traversal time.
   * @param evidence Number of the las iteration of this evidence as counted by AutoTuner.
   */
  virtual void addEvidence(const Configuration &configuration, const Evidence &evidence){};

  /**
   * Selects the next configuration to test or the optimum.
   *
   * A bool is returned indicating whether more tuning steps are required (=true) or the optimum was found (=false).
   * The new configuration can be obtained by getCurrentConfiguration. It is the configuration which is either the next
   * configuration to test (=true) or the optimum (=false).
   *
   * @param configQueue Queue of configurations to be tested. The tuning strategy should edit this queue.
   * @param evidenceCollection All collected evidence until now.
   */
  virtual void optimizeSuggestions(std::vector<Configuration> &configQueue,
                                   const EvidenceCollection &evidenceCollection) = 0;

  /**
   * Reset all internal parameters to the beginning of a new tuning phase.
   * @param iteration Gives the current iteration to the tuning strategy.
   * @param configQueue Queue of configurations to be tested. The tuning strategy should edit this queue.
   * @param evidenceCollection All collected evidence until now.
   */
  virtual void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                     const autopas::EvidenceCollection &evidenceCollection) = 0;

  /**
   * Returns whether this tuning strategy wants to get a LiveInfo object passed before a new tuning phase.
   * @return True, if this tuning strategy wants a LiveInfo object passed before a new tuning phase via
   * receiveLiveInfo().
   */
  [[nodiscard]] virtual bool needsLiveInfo() const { return false; }

  /**
   * Virtual method that subclasses can override to receive the LiveInfo object before a tuning phase if they return
   * true in needsLiveInfo().
   * @param info A new LiveInfo object that has already gathered its information.
   */
  virtual void receiveLiveInfo(const LiveInfo &info){};

  /**
   * Notify the strategy about a configuration that will never be valid and thus can be dropped from any
   * internal storage.
   * @param configuration
   */
  virtual void rejectConfigurationIndefinitely(const Configuration &configuration){};

  /**
   * Indicate whether the strategy needs smoothed values of homogeneity and max density
   * @return
   */
  virtual bool smoothedHomogeneityAndMaxDensityNeeded() const { return false; }
};
}  // namespace autopas
