/**
 * @file TuningTriggerFactoryInfo.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

namespace autopas {

/**
 * Helper struct encapsulating most information needed to build TuningTriggers by the TuningTriggersFactory.
 * This way if only a specific option should be built, not all arguments have to be specified explicitly.
 */
struct TuningTriggerFactoryInfo {
  // Used by multiple triggers, therefore the specific meaning of these variables depends on the trigger that uses them.
  /**
   * Factor on which to trigger when comparing two values derived from iteration samples.
   */
  float triggerFactor;
  /**
   * Number of samples to base trigger decision on.
   */
  unsigned nSamples;
};
}  // namespace autopas
