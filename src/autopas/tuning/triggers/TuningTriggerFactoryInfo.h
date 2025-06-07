/**
 * @file TuningTriggerFactoryInfo.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include <string>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas {

/**
 * Helper struct encapsulating most information needed to build TuningTriggers by the TuningTriggersFactory.
 * This way if only a specific option should be built, not all arguments have to be specified explicitly.
 */
struct TuningTriggerFactoryInfo {
  // Used by multiple triggers
  /**
   * Factor on which to trigger when comparing two values derived from iteration runtimes.
   */
  float factor;
};
}  // namespace autopas
