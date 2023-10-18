/**
 * @file TunerSynchronizer.cpp
 * @author muehlhaeusser
 * @date 09.09.23
 */

#include "TunerSynchronizer.h"

autopas::TunerSynchronizer::TunerSynchronizer(std::set<InteractionTypeOption::Value> &usedInteractionTypes) {
  for (auto &interactionT : usedInteractionTypes) {
    _tuningStates[interactionT] = false;
  }
}

void autopas::TunerSynchronizer::recordTuningState(InteractionTypeOption::Value interactionType, bool stillTuning) {
  _tuningStates[interactionType] = stillTuning;
}

bool autopas::TunerSynchronizer::checkTuningState(autopas::InteractionTypeOption::Value interactionType) {
  if (not _tuningStates[interactionType]) {
    for (const auto &[interactionT, isTuning] : _tuningStates) {
      if (interactionT != interactionType and isTuning) {
        return true;
      }
    }
  }
  return false;
}
