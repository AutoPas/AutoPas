/**
 * @file TuningManager.h
 * @author muehlhaeusser
 * @date 06.02.2026
 */

#pragma once

#include <limits>
#include <memory>
#include <set>
#include <unordered_map>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/tuning/AutoTuner.h"

namespace autopas {

/**
 * The TuningManager is responsible for coordinating the AutoTuner(s)
 * In the case of multiple AutoTuners it ensures that configurations are selected with the same container for all
 * tuners. Otherwise, neighbor lists would need to be rebuilt every time the interaction type (e.g. pairwise and
 * triwise) is switched.
 */
class TuningManager {
 public:
  /**
   * Constructor of the TuningManager.
   * @param autoTunerInfo to get the tuning interval
   */
  explicit TuningManager(const AutoTunerInfo &autoTunerInfo);

  /**
   * Add an AutoTuner to the TuningManager, which will take over ownership.
   * @param interactionType
   * @param tuner A unique_ptr to the new tuner.
   */
  void addAutoTuner(std::unique_ptr<AutoTuner> tuner, InteractionTypeOption::Value interactionType);

  /**
   *
   * @param currentIteration Current LogicHandler iteration number.
   * @param info LiveInfo object to be passed to AutoTuners if relevant.
   * @return true, if at least one AutoTuner is still in tuning state.
   */
  bool tune(size_t currentIteration, const LiveInfo &info);

  /**
   * Add performance metrics for the currently active configuration.
   * @param sampleRebuild Measurement for the rebuild step.
   * @param sampleTraverseParticles Measurement for the traversal step.
   * @param neighborListRebuilt True, if this was a rebuild iteration.
   * @param iteration The current iteration number.
   * @param interactionType
   */
  void addMeasurement(long sampleRebuild, long sampleTraverseParticles, bool neighborListRebuilt, size_t iteration,
                      InteractionTypeOption::Value interactionType) const;

  /**
   * Pass the given rebuild frequency over to all AutoTuners
   * @param rebuildFrequency
   */
  void setRebuildFrequency(double rebuildFrequency);

  /**
   * Reject the current configuration for the interactionType by either skipping over it (indefinitely = false) or
   * removing it entirely (indefinitely = true). Also refreshes the respective AutoTuner's queue and returns the next
   * configuration in line.
   * @param rejectedConfig The configuration to skip/remove
   * @param indefinitely True, if the configuration should also be removed from the search space and with that upcoming
   * tuning phases.
   * @param interactionType Pairwise, triwise, ...
   * @return The next configuration from the queue.
   */
  Configuration rejectConfiguration(const Configuration &rejectedConfig, bool indefinitely,
                                    InteractionTypeOption::Value interactionType);

  /**
   * Force all AutoTuners to start a new tuning phase immediately.
   */
  void forceRetune();

  /**
   * Log the tuning result for a specific interaction type.
   * @param tuningTime Measured time it took for tuning in nanoseconds.
   * @param currentIteration Current LogicHandler iteration number.
   * @param interactionType
   */
  void logTuningResult(long tuningTime, size_t currentIteration, InteractionTypeOption::Value interactionType) const;

  /**
   * Returns whether a rebuild is required during the given iteration, based on AutoTuner(s) providing a new
   * configuration.
   * @param currentIteration Current LogicHandler iteration number.
   * @return True, if a rebuild is necessary due to a configuration change.
   */
  bool requiresRebuilding(size_t currentIteration);

  /**
   *
   * @param currentIteration Current LogicHandler iteration number.
   * @return true, if at least one AutoTuner needs live information for tuning.
   */
  bool needsLiveInfo(size_t currentIteration);

  /**
   * @return True, if the search space of all AutoTuners consists of only one valid configuration.
   */
  bool allSearchSpacesAreTrivial() const;

  /**
   * @return True, if the tuning phase just finished for all AutoTuners.
   */
  bool tuningPhaseJustFinished() const;

  /**
   * Calculate whether we are at the start of a new tuning phase based on the current iteration and tuning interval.
   * @param currentIteration Current LogicHandler iteration number.
   * @return True, if this is the first tuning iteration.
   */
  bool isStartOfTuningPhase(size_t currentIteration) const;

  /**
   * Indicates whether the tuner is in the iteration corresponding to the last sample of the first configuration in the
   * current tuning phase.
   * @param currentIteration Current LogicHandler iteration number.
   * @return True, if this is the last sampling iteration for the first tested config.
   */
  bool inFirstConfigurationLastSample(size_t currentIteration) const;

  /**
   * Get the current configuration for a specific interaction type.
   * @param interactionType
   * @return Configuration
   */
  const Configuration &getCurrentConfig(InteractionTypeOption::Value interactionType) const;

  /**
   * Gets the tuning metric of the AutoTuner of interactionType. (Metric as in time, energy, ...)
   * @param interactionType Get info from the pairwise/triwise AutoTuner.
   * @return TuningMetricOption of the respective AutoTuner.
   */
  const TuningMetricOption &getTuningMetric(InteractionTypeOption::Value interactionType) const;

  /**
   * @return A reference to the map of AutoTuners.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> &getAutoTuners() { return _autoTuners; }

 private:
  /**
   * Let all AutoTuners tune their configuration for the time step ahead.
   */
  void tuneConfigurations(size_t currentIteration);

  /**
   * To be called at the end of a tuning phase.
   * Compares evidence of all active AutoTuners and sets the best configurations leading to an overall optimal time step
   * performance.
   */
  void setOptimalConfigurations();

  /**
   * Find all container options that are part of all currently managed AutoTuner instances.
   */
  std::set<ContainerOption> setCommonContainerOption();

  /**
   * Calculate whether we are about to start a new tuning phase based on the current iteration and tuning interval.
   * @param currentIteration Current LogicHandler iteration number.
   * @return True, if we are within 10 iterations of the next tuning phase.
   */
  bool tuningPhaseAboutToBegin(size_t currentIteration) const;

  /**
   * All AutoTuners used in this instance of AutoPas.
   * There can be up to one per interaction type.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<AutoTuner>> _autoTuners;

  /**
   * Tuning Phase counter.
   */
  size_t _tuningPhase = 0;

  /**
   * New tuning phase starting at multiples of _tuningInterval.
   */
  size_t _tuningInterval = 0;

  /**
   * Number of samples taken for each configuration.
   */
  size_t _maxSamples = 0;

  /**
   * Set to true when the tuning phase just finished for all AutoTuners.
   */
  bool _tuningFinished = true;

  /**
   * Stores the iteration number of the last tuning iteration.
   */
  size_t _lastTuningIteration{std::numeric_limits<size_t>::max()};

  /**
   * Intermediate value to signal that tuning just finished, and we need to select our optimal set-up.
   */
  bool _transitionToOptimalConfigurations = false;

  /**
   * Flag to indicate that a retune has been forced and a new tuning phase should start
   * in the next iteration.
   */
  bool _forceRetunePending{false};
};

}  // namespace autopas