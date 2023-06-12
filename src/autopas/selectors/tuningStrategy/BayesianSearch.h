/**
 * @file BayesianSearch.h
 * @author Jan Nguyen
 * @date 12.06.2019
 */

#pragma once

#include <limits>
#include <map>
#include <set>
#include <unordered_set>

#include "GaussianModel/GaussianProcess.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/selectors/FeatureVectorEncoder.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given
 * feature next.
 */
class BayesianSearch : public TuningStrategyInterface {
  /**
   * The maximum number of attempts to sample an optimum.
   */
  constexpr static size_t maxAttempts = 10;
  /**
   *
   */
  constexpr static double secondsPerMicroseconds = 1. / 1000000.;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumLHSamples number of samples used for prediction while tuning.
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  explicit BayesianSearch(
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      size_t predNumLHSamples = 1000, unsigned long seed = std::random_device()());

  const Configuration &getCurrentConfiguration() const override;

  void removeN3Option(Newton3Option badNewton3Option) override;

  void addEvidence(long time, size_t iteration) override;

  long getEvidence(Configuration configuration) const override;

  void reset(size_t iteration) override;

  bool tune(bool currentInvalid = false) override;

  std::set<ContainerOption> getAllowedContainerOptions() const override;

  bool searchSpaceIsTrivial() const override;

  bool searchSpaceIsEmpty() const override;

  bool smoothedHomogeneityAndMaxDensityNeeded() const override;

 private:
  /**
   * Generate n samples and sort them depending on
   * given acquistion function.
   * @param n numSamples
   * @param af acquisition function
   */
  void sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
  FeatureVectorEncoder _encoder;

  FeatureVector _currentConfig;
  std::vector<FeatureVector> _currentSamples;
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;
  /**
   * Explicitly store traversal times for getEvidence().
   * Refrain from reading the data from GaussianProcesses to maintain abstraction.
   */
  std::unordered_map<Configuration, long, ConfigHash> _traversalTimes;

  Random _rng;
  GaussianProcess _gaussianProcess;
  size_t _maxEvidence;
  AcquisitionFunctionOption _predAcqFunction;
  size_t _predNumLHSamples;
};

}  // namespace autopas
