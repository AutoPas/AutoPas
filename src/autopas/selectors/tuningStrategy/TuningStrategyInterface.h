/**
 * @file BaseTuningStrategy.h
 * @author F. Gratl
 * @date 5/29/19
 */

#pragma once

#include "autopas/selectors/Configuration.h"

namespace autopas {

class TuningStrategyInterface {
  /**
   * Store empirically collected information for a given configuration.
   * @param configuration
   * @param time
   */
  virtual void addEvidence(Configuration configuration, long time) = 0;

  /**
   * Selects the next configuration to test or the optimum.
   *
   * The returned pair first contains a bool that indicates whether more tuning steps are required (=true) or the
   * optimum was found (=false). Second is the configuration which is either the next configuration to test (=true) or
   * the optimum (=false).
   *
   * @return false iff returned configuration is the selected optimum.
   */
  virtual std::pair<bool, Configuration> tune() = 0;
};
}  // namespace autopas