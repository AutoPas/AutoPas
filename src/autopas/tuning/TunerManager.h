/**
 * @file TunerManager.h
 * @author muehlhaeusser
 * @date 06.02.2026
 */

#pragma once

#include <memory>
#include <unordered_map>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/tuning/AutoTuner.h"

namespace autopas {

class TunerManager {
 public:
  TunerManager() = default;

  /**
   * Add an AutoTuner to the TunerManager, which will take over ownership.
   * @param interactionType
   * @param tuner A unique_ptr to the new tuner.
   */
  void addAutoTuner(std::unique_ptr<AutoTuner> tuner, InteractionTypeOption::Value interactionType);

  /**
   * @return A reference to the map of AutoTuners.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> &getAutoTuners() { return _autoTuners; }

 private:
  /**
   * Find all container options that are part of all currently managed AutoTuner instances.
   */
  void setCommonContainerOption();

  /**
   * Constrain all AutoTuners to the given container options and refresh their config queues.
   * @param containerOption
   */
  void applyContainerConstraint(ContainerOption containerOption);

  /**
   * Vector of all allowed container options with configurations for all interaction types.
   */
  std::vector<ContainerOption::Value> _commonContainerOptions;

  /**
   * Index to track the active container option during tuning.
   */
  size_t _currentContainerIndex = 0;

  /**
   * All AutoTuners used in this instance of AutoPas.
   * There can be up to one per interaction type.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> _autoTuners;
};

}  // namespace autopas