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

/**
 * The TunerManager is responsible for coordinating the AutoTuner(s)
 * In the case of multiple AutoTuners it ensures that configurations are selected with the same container for all
 * tuners. Otherwise, neighbor lists would need to be rebuilt every time the interaction type (e.g. pairwise and
 * triwise) is switched.
 */
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
   * @param currentIteration
   */
  void updateAutoTuners(size_t currentIteration);

  bool tune(size_t currentIteration, const LiveInfo &info);

  bool needsLiveInfo(size_t currentIteration);

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
   * Force all AutoTuners to start a new tuning phase immediately.
   */
  void forceRetune();

  Configuration rejectConfiguration(const Configuration &rejectedConfig, bool indefinitely, InteractionTypeOption::Value interactionType);

  /**
   * @return A reference to the map of AutoTuners.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> &getAutoTuners() { return _autoTuners; }

 private:
  /**
   * Find all container options that are part of all currently managed AutoTuner instances.
   */
  std::set<ContainerOption> setCommonContainerOption();

  /**
   * Constrain all AutoTuners to the given container options and refresh their config queues.
   * @param containerOption
   */
  void applyContainerConstraint(std::optional<ContainerOption> containerOption);

  /**
   * Let all AutoTuners tune their configuration for the time step ahead.
   */
  void tuneConfigurations();

  void setOptimalConfigurations();

  /**
   * All AutoTuners used in this instance of AutoPas.
   * There can be up to one per interaction type.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> _autoTuners;

  /**
   * New tuning phase starting at multiples of _tuningInterval.
   */
  size_t _tuningInterval = 0;

  /**
   * Set to true when the tuning phase just finished for all AutoTuners.
   */
  bool _tuningFinished = true;

  size_t _lastTuningIteration{std::numeric_limits<size_t>::max()};

  /**
   * Initialize
   */
  bool _tuningJustFinished = false;
};

}  // namespace autopas