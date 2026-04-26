/**
 * @file ReinforcementLearning.h
 * @author P. Metscher
 * @date 07.05.2025
 */

#pragma once

#include <unordered_map>

#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {
/**
 * Using reinforcement learning to determine the best configuration.
 *
 * This class is a reimplementation based on L. Laumeyer's implementation of the reinforcement learning algorithm
 * selection.
 */
class ReinforcementLearning final : public TuningStrategyInterface {
 public:
  /**
   * Create a reinforcement learning tuning strategy.
   * @param searchSpace The search space to use.
   * @param learningRate The learning rate for the reinforcement learning algorithm.
   * @param discountFactor The discount factor for the reinforcement learning algorithm.
   * @param numRandomExplorations Optional: The number of random explorations to perform.
   */
  explicit ReinforcementLearning(const std::set<Configuration> &searchSpace, const double learningRate,
                                 const double discountFactor, const size_t numRandomExplorations = 5);

  ~ReinforcementLearning() override = default;

  /**
   * Get this object's associated TuningStrategyOption type.
   * @return TuningStrategyOption::reinforcementLearning
   */
  TuningStrategyOption getOptionType() const override;

  /**
   * Notifies the strategy about empirically collected information for the given configuration.
   *
   * This function is required for updating the state of the reinforcement learning algorithm.
   *
   * @param configuration Configuration used to obtain the evidence.
   * @param evidence Measurement and when it was taken.
   */
  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  /**
   * Optimizes the queue of configurations to process.
   *
   * This function is called once before each iteration in a tuning phase so all tuning strategies can give their
   * input on which configuration to try next. This is done by reordering configQueue so that the next configuration
   * to try is at the end (FIFO).
   *
   * @param configQueue Queue of configurations to be tested. The tuning strategy should edit this queue.
   * @param evidenceCollection All collected evidence until now.
   * @return boolean value to signal if the tuning strategy has intentionally wiped the config queue
   */
  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  /**
   * Reset all internal parameters to the beginning of a new tuning phase.
   *
   * This can also mean to reorder the configQueue to some initially expected state.
   *
   * @param iteration Gives the current iteration to the tuning strategy.
   * @param tuningPhase Gives the current tuning phase to the tuning strategy.
   * @param configQueue Queue of configurations to be tested. The tuning strategy should edit this queue.
   * @param evidenceCollection All collected evidence until now.
   * @return boolean value to signal if the tuning strategy has intentionally wiped the config queue
   */
  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const EvidenceCollection &evidenceCollection) override;

 private:
  /**
   * Store the learning rate.
   */
  double _learningRate;

  /**
   * Store the discount factor.
   */
  double _discountFactor;

  /**
   * Store the number of random explorations.
   *
   * This defines how many random explorations will take place every tuning phase.
   */
  size_t _numRandomExplorations = 5;

  /**
   * Store the state of the reinforcement learning algorithm.
   *
   * This is a map of configurations and double weights.
   */
  std::unordered_map<Configuration, double, ConfigHash> _state;

  /**
   * Store the best evidence in this tuning phase.
   *
   * @note This is a long and not a double as it is only directly compared to evidence.value.
   */
  long _bestEvidence = std::numeric_limits<long>::max();

  /**
   * Store the best configuration during this tuning phase.
   */
  Configuration _bestConfiguration;

  /**
   * Store if this tuning phase is the first tuning phase.
   */
  bool _firstTuningPhase = true;

  /**
   * Store if the program is in an exploration phase.
   */
  bool _explorationPhase = true;
};
}  // namespace autopas