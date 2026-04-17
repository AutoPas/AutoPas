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
  bool requiresRebuilding(size_t currentIteration);

  /**
   * @return True, if the search space of all AutoTuners consists of only one valid configuration.
   */
  bool allSearchSpacesAreTrivial() const;

  /**
   * @return True, if the tuning phase just finished for all AutoTuners.
   */
  bool tuningPhaseJustFinished() const;

  /**
   * Log the tuning result for a specific interaction type.
   */
  void logTuningResult(long tuningTime, size_t currentIteration, InteractionTypeOption::Value interactionType) const;

  /**
   * Add a performance measurement to the specified tuner.
   */
  void addMeasurement(long sampleRebuild, long sampleTraverseParticles, bool neighborListRebuilt, size_t iteration,
                      InteractionTypeOption::Value interactionType);

  const TuningMetricOption &getTuningMetric(InteractionTypeOption::Value interactionType) const;

  /**
   * Get the current configuration for a specific interaction type.
   */
  const Configuration &getCurrentConfig(InteractionTypeOption::Value interactionType) const;

  /**
   * Force all AutoTuners to start a new tuning phase immediately.
   */
  void forceRetune();

  Configuration rejectConfiguration(const Configuration &rejectedConfig, bool indefinitely, size_t currentIteration,
                                    InteractionTypeOption::Value interactionType);

  /**
   * @return A reference to the map of AutoTuners.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> &getAutoTuners() { return _autoTuners; }

  bool inFirstTuningIteration(size_t currentIteration) const;

 private:
  /**
   * Find all container options that are part of all currently managed AutoTuner instances.
   */
  std::set<ContainerOption> setCommonContainerOption();

  /**
   * Let all AutoTuners tune their configuration for the time step ahead.
   */
  void tuneConfigurations(size_t currentIteration);

  void setOptimalConfigurations(size_t currentIteration);

  bool tuningPhaseAboutToBegin(size_t currentIteration) const;

  bool isStartOfTuningPhase(size_t currentIteration) const;
  /**
   * All AutoTuners used in this instance of AutoPas.
   * There can be up to one per interaction type.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> _autoTuners;

  size_t _tuningPhase = 0;
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
  bool _transitionToOptimalConfigurations = false;

  /**
   * Flag to indicate that a retune has been forced and a new tuning phase should start
   * in the next iteration.
   */
  bool _forceRetunePending{false};
};

}  // namespace autopas