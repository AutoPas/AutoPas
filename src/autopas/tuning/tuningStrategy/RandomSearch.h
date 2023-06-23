/**
 * @file RandomSearch.h
 * @author Jan Nguyen
 * @date 10.07.2019
 */

#pragma once

#include <map>
#include <set>

#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Randomly test a given number of configurations and select the fastest.
 */
class RandomSearch : public TuningStrategyInterface {
 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  explicit RandomSearch(
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1. / 3., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), size_t maxEvidence = 10,
      unsigned long seed = std::random_device()());

  const Configuration &getCurrentConfiguration() const override;

  void removeN3Option(Newton3Option badNewton3Option) override;

  void addEvidence(long time, size_t iteration) override;

  long getEvidence(Configuration configuration) const override;

  void reset(size_t iteration) override;

  bool tune(bool currentInvalid = false) override;

  std::set<ContainerOption> getAllowedContainerOptions() const override;

  bool searchSpaceIsTrivial() const override;

  bool searchSpaceIsEmpty() const override;

  bool smoothedHomogeneityAndMaxDensityNeeded() const override { return false; }

 private:
  void selectOptimalConfiguration();

  std::set<ContainerOption> _containerOptions;
  std::set<TraversalOption> _traversalOptions;
  std::set<LoadEstimatorOption> _loadEstimatorOptions;
  std::set<DataLayoutOption> _dataLayoutOptions;
  std::set<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;

  Configuration _currentConfig;
  std::unordered_map<Configuration, size_t, ConfigHash> _traversalTimes;
  std::set<Configuration> _invalidConfigs;

  Random _rng;
  size_t _maxEvidence;
  size_t _searchSpaceSizeNoCSF;
  size_t _numTestedConfigsNoCSF;
};

}  // namespace autopas
