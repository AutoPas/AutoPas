/**
 * @file TuningTriggerFactory.cpp
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#include "TuningTriggerFactory.h"

#include "autopas/options/TuningTriggerOption.h"
#include "autopas/tuning/triggers/StaticSimpleTrigger.h"
#include "autopas/tuning/triggers/TimeBasedAverageTrigger.h"
#include "autopas/tuning/triggers/TimeBasedRegressionTrigger.h"
#include "autopas/tuning/triggers/TimeBasedSimpleTrigger.h"
#include "autopas/tuning/triggers/TimeBasedSplitTrigger.h"

namespace autopas::TuningTriggerFactory {

std::unique_ptr<TuningTriggerInterface> generateTuningTrigger(TuningTriggerOption tuningTriggerOption,
                                                              const TuningTriggerFactoryInfo &info) {
  std::unique_ptr<TuningTriggerInterface> tuningTrigger = nullptr;
  switch (static_cast<TuningTriggerOption>(tuningTriggerOption)) {
    case TuningTriggerOption::staticSimple: {
      tuningTrigger = std::make_unique<StaticSimpleTrigger>();
      break;
    }
    case TuningTriggerOption::timeBasedSimple: {
      tuningTrigger = std::make_unique<TimeBasedSimpleTrigger>(info.triggerFactor);
      break;
    }
    case TuningTriggerOption::timeBasedAverage: {
      tuningTrigger = std::make_unique<TimeBasedAverageTrigger>(info.triggerFactor, info.nSamples);
      break;
    }
    case TuningTriggerOption::timeBasedSplit: {
      tuningTrigger = std::make_unique<TimeBasedSplitTrigger>(info.triggerFactor, info.nSamples);
      break;
    }
    case TuningTriggerOption::timeBasedRegression: {
      tuningTrigger = std::make_unique<TimeBasedRegressionTrigger>(info.triggerFactor, info.nSamples);
      break;
    }
    default: {
      utils::ExceptionHandler::exception("AutoPas::generateTuningTrigger: Unknown tuning trigger {}!",
                                         tuningTriggerOption);
      break;
    }
  }
  return tuningTrigger;
}
}  // namespace autopas::TuningTriggerFactory