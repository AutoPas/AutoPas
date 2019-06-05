/**
 * @file TuningStrategyInterface.h
 * @author F. Gratl
 * @date 5/29/19
 */

#pragma once

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
   */
  virtual void addEvidence(long time) = 0;

  /**
   * Returns the currently selected configuration object.
   * @return
   */
  virtual Configuration getCurrentConfiguration() = 0;

  /**
   * Selects the next configuration to test or the optimum.
   *
   * A bool is returned indicating whether more tuning steps are required (=true) or the optimum was found (=false).
   * The new configuration can be obtained by getCurrentConfiguration. It is the configuration which is either the next
   * configuration to test (=true) or the optimum (=false).
   *
   * @return false iff new configuration is the selected optimum.
   */
  virtual bool tune() = 0;

  /**
   * Reset all internal parameters to the beginning of a new tuning cycle.
   */
  virtual void reset() = 0;

  /**
   * Returns all container options the strategy might choose.
   * @return
   */
  virtual std::set<ContainerOption> getAllowedContainerOptions() = 0;

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
  virtual bool searchSpaceOneOption() = 0;

  /**
   * Indicate whether the search space collapsed to be empty.
   * @return
   */
  virtual bool searchSpaceEmpty() = 0;
};
}  // namespace autopas