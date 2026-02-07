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
  /**
   * Constructor of the TunerManager.
   * @param autoTunerInfo to get the tuning interval
   */
  explicit TunerManager(const AutoTunerInfo &autoTunerInfo);

  /**
   * Add an AutoTuner to the TunerManager, which will take over ownership.
   * @param interactionType
   * @param tuner A unique_ptr to the new tuner.
   */
  void addAutoTuner(std::unique_ptr<AutoTuner> tuner, InteractionTypeOption::Value interactionType);

  /**
   * Prepare AutoTuners for the next time step. Bumps internal counters and tunes configurations.
   */
  void updateAutoTuners();

  /**
   * Reject given configuration from the AutoTuners config queue and search space (if indefinitely==true). Also handles
   * the case if the rejected config was the last one for the currently selected container option.
   * @param configuration
   * @param indefinitely
   * @return
   */
  Configuration rejectConfig(const Configuration &configuration, bool indefinitely);

  /**
   * @return True, if a rebuild is necessary due to a configuration change.
   */
  bool requiresRebuilding();

  /**
   * @return True, if the search space of all AutoTuners consists of only one valid configuration.
   */
  bool allSearchSpacesAreTrivial() const;

  /**
   * @return True, if the tuning phase just finished for all AutoTuners.
   */
  bool tuningPhaseJustFinished() const;

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
   * Increment iteration counters also for all AutoTuners. Should be called exactly once per time step.
   */
  void bumpTunerCounters();

  /**
   * Let all AutoTuners tune their configuration for the time step ahead.
   */
  void tuneConfigurations();

  /**
   * Store the best result for the current container option.
   */
  void captureCurrentContainerPerformance();

  /**
   * Find the container option with the best results in _containerResults.
   * Apply the chosen container option to all AutoTuners and tune once to push the best configuration to the config
   * queue.
   */
  void selectBestContainer();

  /**
   * Marks the current container as invalid as at least one autotuner couldn't find a valid configuration with it.
   * @return true if a new container is available, false if we ran out of options.
   */
  void rejectCurrentContainer();

  /**
   * All AutoTuners used in this instance of AutoPas.
   * There can be up to one per interaction type.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> _autoTuners;

  /**
   * Vector of all allowed container options with configurations for all interaction types.
   */
  std::vector<ContainerOption::Value> _commonContainerOptions;

  /**
   * Store the tuning results for each container option of the last tuning phase.
   */
  std::unordered_map<ContainerOption::Value, long> _containerResults;

  /**
   * Index to track the active container option during tuning.
   */
  size_t _currentContainerIndex = 0;

  /**
   * Iteration counter. Starts at max value, so ++_iteration == 0 in the first time step.
   */
  size_t _iteration = std::numeric_limits<size_t>::max();

  /**
   * New tuning phase starting at multiples of _tuningInterval.
   */
  size_t _tuningInterval = 0;

  /**
   * Set to true when the tuning phase just finished for all AutoTuners.
   */
  bool _tuningJustFinished = false;
};

}  // namespace autopas