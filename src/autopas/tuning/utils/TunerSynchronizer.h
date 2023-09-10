/**
 * @file TunerSynchronizer.h
 * @author muehlhaeusser
 * @date 09.09.23
 */

#pragma once

#include <set>
#include "autopas/options/InteractionTypeOption.h"

namespace autopas {
/**
 * Track synchronization states of multiple autotuners.
 * To be used such that new tuning phases start at the same time for all autotuners.
 */
class TunerSynchronizer {
 public:


  TunerSynchronizer() = default;

  /**
   * Constructor
   */
  explicit TunerSynchronizer(std::set<InteractionTypeOption::Value> &usedInteractionTypes);

  /**
   * Add a new interaction type to keep track of.
   * @param interactionType interaction type to add for synchronization
   * @param isTuning initial tuning status (default to true)
   */
  void addInteractionType(InteractionTypeOption::Value interactionType, const bool isTuning = true) {_tuningStates[interactionType] = isTuning;}

  /**
   * To be called after configuration was selected and tuning status determined.
   * @param interactionType
   * @param stillTuning
   */
  void recordTuningState(InteractionTypeOption::Value interactionType, bool stillTuning);

  /**
   *
   * @param interactionType
   * @return bool, whether other tuners are still tuning
   */
  bool checkState(InteractionTypeOption::Value interactionType);

 private:
  /**
   * Dictionary of interaction type and tuning status of the corresponding tuner.
   */
  std::unordered_map<InteractionTypeOption::Value, bool> _tuningStates{};
};

}  // namespace autopas