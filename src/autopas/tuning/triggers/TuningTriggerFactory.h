/**
 * @file TuningTriggerFactory.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include "TuningTriggerFactoryInfo.h"
#include "autopas/options/TuningTriggerOption.h"
#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas::TuningTriggerFactory {
/**
 * Generates a new Tuning Trigger object.
 * @param tuningTriggerOption
 * @param info
 * @return Pointer to the tuning trigger object, or null pointer if an exception was suppressed.
 */
std::unique_ptr<TuningTriggerInterface> generateTuningTrigger(TuningTriggerOption tuningTriggerOption,
                                                              const TuningTriggerFactoryInfo &info);
}  // namespace autopas::TuningTriggerFactory
