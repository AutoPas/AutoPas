/**
 * @file BayesianClusterSearch.h
 * @author Jan Nguyen
 * @date 12.04.20
 */

#pragma once

#include <limits>
#include <map>
#include <set>
#include <unordered_set>

#include "FullSearch.h"
#include "GaussianModel/GaussianCluster.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/tuning/utils/FeatureVectorEncoder.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time while fixing discrete variables corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given feature next.
 */
class BayesianClusterSearch : public TuningStrategyInterface {
  /**
   * The maximum number of attempts to sample an optimum.
   */
  constexpr static size_t maxAttempts = 10;
  /**
   * The factor for conversion from seconds to microseconds.
   */
  constexpr static double secondsPerMicroseconds = 1. / 1000000.;
  /**
   * Dimensions of the continuous tuples.
   * FeatureVector continuous dimensions + 1 (time)
   */
  constexpr static size_t continuousDims = FeatureVectorEncoder::tunableContinuousDims + 1;
  /**
   * Fixed noise
   */
  constexpr static double sigma = 0.01;

  /**
   * Iteration numbers get scaled down by a multiple of max evidence to ensure not all kernels become 0.
   */
  constexpr static double iterationScalePerMaxEvidence = 5.;

  /**
   * Distance between evidence is expected not to exceed a upper threshold.
   * We chose the threshold 7 because: exp(-r) < 0.1% for r > 7 but this
   * can be increased if necessary.
   */
  constexpr static double suggestedMaxDistance = 7.;

 public:
  /**
   * Constructor
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param maxEvidence Stop tuning after given number of evidence provided.
   * @param predAcqFunction Acquisition function used for prediction while tuning.
   * @param outputSuffix Suffix for output logger.
   * @param predNumLHSamples Number of latin-hypercube-samples used to find a evidence with high predicted acquisition
   * @param seed Seed of random number generator (should only be used for tests)
   */
  explicit BayesianClusterSearch(
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      const std::string &outputSuffix = "", size_t predNumLHSamples = 50, unsigned long seed = std::random_device()());

  ~BayesianClusterSearch() override;

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
  /**
   * Generate n samples and predict their corresponding
   * acquisition function. The result is sorted depending on
   * the given acquisition function.
   * @param n numSamples
   * @param af acquisition function
   */
  void sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  /**
   * If allowed options are changed this functions should be called
   * to update encoder and clusters.
   */
  void updateOptions();

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
  FeatureVectorEncoder _encoder;

  FeatureVector _currentConfig;
  /**
   * Currently sampled vectors and corresponding acquisition values.
   */
  std::vector<GaussianModelTypes::VectorPairDiscreteContinuous> _currentAcquisitions;
  /**
   * Configurations marked invalid.
   */
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;
  /**
   * Explicitly store traversal times for getEvidence().
   * Refrain from reading the data from GaussianProcesses to maintain abstraction.
   */
  std::unordered_map<Configuration, long, ConfigHash> _traversalTimes;
  /**
   * Random engine.
   */
  Random _rng;
  /**
   * Stochastic model used for predictions.
   */
  GaussianCluster _gaussianCluster;
  /**
   * Function to generate neighbours with weight of given vector.
   */
  GaussianModelTypes::NeighbourFunction _neighbourFun;
  /**
   * Maximum number of evidences to collect after which an optimum is declared.
   */
  const size_t _maxEvidence;
  /**
   * Acquisition function used to predict informational gain.
   */
  const AcquisitionFunctionOption _predAcqFunction;
  /**
   * Number of latin-hypercube-samples used to find a evidence with high predicted acquisition.
   */
  const size_t _predNumLHSamples;
  /**
   * Flag indicating that we are in the first tuning phase.
   */
  bool _firstTuningPhase;
  /**
   * Iteration of last added evidence or reset.
   */
  size_t _currentIteration;
  /**
   * Iteration numbers get scaled down by a multiple of max evidence to ensure not all kernels become 0.
   */
  const double _iterationScale;
  /**
   * Number of evidence provided in current tuning phase.
   */
  size_t _currentNumEvidence;
  /**
   * Lowest time of evidence in current tuning phase.
   */
  long _currentOptimalTime;
  /**
   * Configuration with lowest time in current tuning phase.
   */
  FeatureVector _currentOptimalConfig;

  /**
   * FullSearch used in the first tuning phase.
   */
  FullSearch _fullSearch;
};
}  // namespace autopas
