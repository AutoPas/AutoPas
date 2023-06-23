/**
 * @file SetSearchSpaceBasedTuningStrategy.h
 * @author Julian Pelloth
 * @date 29.04.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Super class for tuning strategies with a set based search space.
 */
class SetSearchSpaceBasedTuningStrategy : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the SetSearchSpaceBasedTuningStrategy that generates the search space from the allowed options.
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   */
  SetSearchSpaceBasedTuningStrategy(const std::set<ContainerOption> &allowedContainerOptions,
                                    const std::set<double> &allowedCellSizeFactors,
                                    const std::set<TraversalOption> &allowedTraversalOptions,
                                    const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                                    const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                                    const std::set<Newton3Option> &allowedNewton3Options);

  /**
   * Constructor for the SetSearchSpaceBasedTuningStrategy that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   */
  explicit SetSearchSpaceBasedTuningStrategy(std::set<Configuration> allowedConfigurations);

  [[nodiscard]] std::set<ContainerOption> getAllowedContainerOptions() const override;

  [[nodiscard]] const Configuration &getCurrentConfiguration() const override;

  [[nodiscard]] bool searchSpaceIsTrivial() const override;

  [[nodiscard]] bool searchSpaceIsEmpty() const override;

  [[nodiscard]] bool smoothedHomogeneityAndMaxDensityNeeded() const override;

  void removeN3Option(Newton3Option badNewton3Option) override;

 protected:
  /**
   * Fills the search space with the cartesian product of the given options (minus invalid combinations).
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   */
  void populateSearchSpace(const std::set<ContainerOption> &allowedContainerOptions,
                           const std::set<double> &allowedCellSizeFactors,
                           const std::set<TraversalOption> &allowedTraversalOptions,
                           const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions,
                           const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                           const std::set<Newton3Option> &allowedNewton3Options);

  /**
   * Finds the optimal configuration in a given search space.
   * @param searchSpace Map<Configuration, runTime> to search.
   * @return Iterator to optimal (=fastest) configuration.
   */
  static std::unordered_map<autopas::Configuration, long, autopas::ConfigHash>::const_iterator getOptimum(
      const std::unordered_map<Configuration, long, ConfigHash> &searchSpace);

  /**
   * Contains every allowed container option.
   */
  std::set<ContainerOption> _allowedContainerOptions{};
  /**
   * Contains every configuration.
   */
  std::set<Configuration> _searchSpace{};
  /**
   * Contains the corrent configuration.
   */
  std::set<Configuration>::iterator _currentConfig{};
};
}  // namespace autopas
