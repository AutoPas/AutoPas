/**
 * @file ReinforcementLearning.h
 * @author L. Laumeyer
 * @date 17.05.2022
 */

#pragma once

#include <set>
#include <sstream>
#include <utility>

#include "SetSearchSpaceBasedTuningStrategy.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 *
 */
class ReinforcementLearning : public SetSearchSpaceBasedTuningStrategy {
 public:
  /**
   * Constructor for the FullSearch that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   */
  ReinforcementLearning(const std::set<ContainerOption> &allowedContainerOptions,
                        const std::set<double> &allowedCellSizeFactors,
                        const std::set<TraversalOption> &allowedTraversalOptions,
                        const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                        const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                        const std::set<Newton3Option> &allowedNewton3Options)
      : SetSearchSpaceBasedTuningStrategy(allowedContainerOptions, allowedCellSizeFactors, allowedTraversalOptions,
                                          allowedLoadEstimatorOptions, allowedDataLayoutOptions,
                                          allowedNewton3Options) {}

  /**
   * Constructor for the ReinforcementLearning that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit ReinforcementLearning(std::set<Configuration> allowedConfigurations)
      : SetSearchSpaceBasedTuningStrategy(std::move(allowedConfigurations)) {}

  inline void addEvidence(long time, size_t iteration) override {
    // TODO
  }

  inline long getEvidence(Configuration configuration) const override { _traversalTimes[*_currentConfig] = time; }

  inline const Configuration &getCurrentConfiguration() const override { return *_currentConfig; }

  inline void reset(size_t iteration) override { _traversalTimes.clear(); }

  inline bool tune(bool = false) override;

 private:
  inline void selectOptimalConfiguration();

  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
};

bool ReinforcementLearning::tune(bool) {
  // TODO
}

}  // namespace autopas
