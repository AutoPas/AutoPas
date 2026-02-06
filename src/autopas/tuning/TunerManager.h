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

  void addAutoTuner(InteractionTypeOption::Value interactionType, std::unique_ptr<AutoTuner> tuner);

 private:
  void setCommonContainerOption();

  void applyContainerConstraint(ContainerOption::Value container);

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