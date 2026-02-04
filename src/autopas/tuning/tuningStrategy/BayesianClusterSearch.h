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

#include "GaussianModel/GaussianCluster.h"
#include "TuningStrategyInterface.h"
#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
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
   * @param interactionType
   * @param allowedContainerOptions
   * @param allowedCellSizeFactors
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedThreadCounts
   * @param maxEvidence Stop tuning after given number of evidence provided.
   * @param predAcqFunction Acquisition function used for prediction while tuning.
   * @param outputSuffix Suffix for output logger.
   * @param predNumLHSamples Number of latin-hypercube-samples used to find a evidence with high predicted acquisition
   * @param seed Seed of random number generator (should only be used for tests)
   */
  explicit BayesianClusterSearch(
      const InteractionTypeOption &interactionType,
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(),
      const NumberSet<int> &allowedThreadCounts = NumberSetFinite<int>({ autopas::Configuration::ThreadCountNoTuning }),
      size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      const std::string &outputSuffix = "", size_t predNumLHSamples = 50, unsigned long seed = std::random_device()());

  ~BayesianClusterSearch() override;

  TuningStrategyOption getOptionType() const override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  void rejectConfiguration(const Configuration &configuration, bool indefinitely) override;

  /**
   * Indicate if the search space is empty.
   * @return
   */
  bool searchSpaceIsEmpty() const;

  bool needsDomainSimilarityStatistics() const override { return false; }

 private:
  /**
   * Generate n samples and predict their corresponding
   * acquisition function. The result is sorted depending on
   * the given acquisition function.
   * @param n numSamples
   * @param af acquisition function
   * @return Vector of acquisitions.
   */
  std::vector<autopas::GaussianModelTypes::VectorPairDiscreteContinuous> sampleAcquisitions(
      size_t n, AcquisitionFunctionOption af);

  /**
   * If allowed options are changed this functions should be called
   * to update encoder and clusters.
   */
  void updateOptions();

  const InteractionTypeOption _interactionType;

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
  std::unique_ptr<NumberSet<int>> _threadCounts;
  FeatureVectorEncoder _encoder;

  /**
   * Currently sampled vectors and corresponding acquisition values.
   */
  std::vector<GaussianModelTypes::VectorPairDiscreteContinuous> _currentAcquisitions;
  /**
   * Configurations marked invalid.
   */
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;
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
};
}  // namespace autopas
