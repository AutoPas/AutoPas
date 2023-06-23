/**
 * @file FullSearch.h
 * @author F. Gratl
 * @date 29.05.2019
 */

#pragma once

#include <set>
#include <sstream>
#include <utility>

#include "SetSearchSpaceBasedTuningStrategy.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/tuning/selectors/OptimumSelector.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Exhaustive full search of the search space by testing every applicable configuration and then selecting the optimum.
 */
class FullSearch : public SetSearchSpaceBasedTuningStrategy {
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
  FullSearch(const std::set<ContainerOption> &allowedContainerOptions, const std::set<double> &allowedCellSizeFactors,
             const std::set<TraversalOption> &allowedTraversalOptions,
             const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
             const std::set<DataLayoutOption> &allowedDataLayoutOptions,
             const std::set<Newton3Option> &allowedNewton3Options);

  /**
   * Constructor for the FullSearch that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit FullSearch(std::set<Configuration> allowedConfigurations);

  void addEvidence(long time, size_t iteration) override;

  long getEvidence(Configuration configuration) const override;

  const Configuration &getCurrentConfiguration() const override;

  void reset(size_t iteration) override;

  bool tune(bool = false) override;

 protected:
  /**
   * Selects the optimal configuration.
   */
  void selectOptimalConfiguration();

  /**
   * Stores the time for each configuration that was already tested.
   */
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
};
}  // namespace autopas
