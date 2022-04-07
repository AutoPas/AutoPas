/**
 * @file TuningStrategyInterface.h
 * @author F. Gratl
 * @date 29.05.2019
 */

#pragma once

#include "LiveInfo.h"
#include "autopas/selectors/Configuration.h"

namespace autopas {

/**
 * Interface for tuning strategies for the auto tuner.
 */
class TuningStrategyInterface {
 public:
  virtual ~TuningStrategyInterface() = default;

  /**
   * Store empirically collected information for the current configuration.
   * @param time Measured traversal time.
   * @param iteration Number of the las iteration of this evidence as counted by AutoTuner.
   */
  virtual void addEvidence(long time, size_t iteration) = 0;

  /**
   * Get the stored time for the given configuration
   * @param configuration
   * @return
   */
  virtual long getEvidence(Configuration configuration) const = 0;

  /**
   * Returns the currently selected configuration object.
   * @return
   */
  virtual const Configuration &getCurrentConfiguration() const = 0;

  /**
   * Selects the next configuration to test or the optimum.
   *
   * A bool is returned indicating whether more tuning steps are required (=true) or the optimum was found (=false).
   * The new configuration can be obtained by getCurrentConfiguration. It is the configuration which is either the next
   * configuration to test (=true) or the optimum (=false).
   *
   * @param currentInvalid Tells the tune() function that the currently selected configuration is invalid. This can be
   * used to avoid getting stuck in an invalid optimum.
   * @return false iff new configuration is the selected optimum.
   */
  virtual bool tune(bool currentInvalid = false) = 0;

  /**
   * Reset all internal parameters to the beginning of a new tuning cycle.
   * @param iteration Gives the current iteration to the tuning strategy.
   */
  virtual void reset(size_t iteration) = 0;

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
  virtual void receiveLiveInfo(const LiveInfo& info){};

  /**
   * Returns all container options the strategy might choose.
   * @return
   */
  virtual std::set<ContainerOption> getAllowedContainerOptions() const = 0;

  /**
   * Removes all configurations with the given newton 3 option from the search space.
   *
   * If the current configuration is removed, it is set to the next not-removed one.
   */
  virtual void removeN3Option(Newton3Option) = 0;

  /**
   * Indicate whether the search space collapsed to only one option.
   * @return
   */
  virtual bool searchSpaceIsTrivial() const = 0;

  /**
   * Indicate whether the search space collapsed to be empty.
   * @return
   */
  virtual bool searchSpaceIsEmpty() const = 0;
};
}  // namespace autopas
