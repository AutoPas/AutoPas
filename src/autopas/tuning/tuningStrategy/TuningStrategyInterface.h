/**
 * @file TuningStrategyInterface.h
 * @author F. Gratl
 * @date 29.05.2019
 */

#pragma once

#include <functional>
#include <optional>

#include "autopas/options/TuningStrategyOption.h"
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
   * Get this object's associated TuningStrategyOption type.
   * @return TuningStrategyOption
   */
  virtual TuningStrategyOption getOptionType() const = 0;

  /**
   * Notifies the strategy about empirically collected information for the given configuration.
   *
   * All evidence is stored centrally in the AutoTuner and its EvidenceCollection is passed to the tuning strategies
   * during optimization.
   *
   * Implementing this function is only necessary if the tuning strategy processes evidence differently
   * than EvidenceCollection.
   *
   * @param configuration Configuration used to obtain the evidence.
   * @param evidence Measurement and when it was taken.
   */
  virtual void addEvidence(const Configuration &configuration, const Evidence &evidence){};

  /**
   * Optimizes the queue of configurations to process.
   *
   * This function is called once before each iteration in a tuning phase so all tuning strategies can give their
   * input on which configuration to try next. This is done by reordering configQueue so that the next configuration
   * to try is at the end (FIFO).
   *
   * @param configQueue Queue of configurations to be tested. The tuning strategy should edit this queue.
   * @param evidenceCollection All collected evidence until now.
   * @param intentionalConfigWipe Optional flag that is used by tuning strategies to signal they have intentionally
   * wiped the config queue
   */
  virtual void optimizeSuggestions(
      std::vector<Configuration> &configQueue, const EvidenceCollection &evidenceCollection,
      std::optional<std::reference_wrapper<bool>> intentionalConfigWipe = std::nullopt) = 0;

  /**
   * Reset all internal parameters to the beginning of a new tuning phase.
   *
   * This can also mean to reorder the configQueue to some initially expected state.
   *
   * @param iteration Gives the current iteration to the tuning strategy.
   * @param tuningPhase Gives the current tuning phase to the tuning strategy.
   * @param configQueue Queue of configurations to be tested. The tuning strategy should edit this queue.
   * @param evidenceCollection All collected evidence until now.
   * @param intentionalConfigWipe Optional flag that is used by tuning strategies to signal they have intentionally
   * wiped the config queue
   */
  virtual void reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                     const autopas::EvidenceCollection &evidenceCollection,
                     std::optional<std::reference_wrapper<bool>> intentionalConfigWipe = std::nullopt) = 0;

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
   * Notify the strategy about a configuration that is (currently) invalid and thus can potentially be dropped from some
   * internal storage.
   * @param configuration
   * @param indefinitely Whether the given configuration will never be valid
   */
  virtual void rejectConfiguration(const Configuration &configuration, bool indefinitely){};

  /**
   * Indicate whether the strategy needs smoothed values of homogeneity and max density
   * @return
   */
  virtual bool needsSmoothedHomogeneityAndMaxDensity() const { return false; }

  /**
   * Method to pass smoothed homogeneity and the maximal density to the tuning strategy.
   * @param homogeneity
   * @param maxDensity
   */
  virtual void receiveSmoothedHomogeneityAndMaxDensity(double homogeneity, double maxDensity){};
};
}  // namespace autopas
