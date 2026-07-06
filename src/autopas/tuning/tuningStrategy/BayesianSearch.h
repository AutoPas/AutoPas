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
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/utils/FeatureVectorEncoder.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Assume that the stochastic distribution of the execution time corresponds
 * to a Gaussian Process. This allows to estimate the 'gain' of testing a given
 * feature next.
 */
class BayesianSearch final : public TuningStrategyInterface {
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
   * @param interactionType
   * @param allowedContainerOptions
   * @param allowedTraversalOptions
   * @param allowedLoadEstimatorOptions
   * @param allowedDataLayoutOptions
   * @param allowedNewton3Options
   * @param allowedCellSizeFactors
   * @param allowedVecPatternOptions
   * @param predAcqFunction acquisition function used for prediction while tuning.
   * @param predNumLHSamples number of samples used for prediction while tuning.
   * @param maxEvidence stop tuning after given number of evidence provided.
   * @param seed seed of random number generator (should only be used for tests)
   */
  explicit BayesianSearch(
      const InteractionTypeOption &interactionType,
      const std::set<ContainerOption> &allowedContainerOptions = ContainerOption::getAllOptions(),
      const NumberSet<double> &allowedCellSizeFactors = NumberInterval<double>(1., 2.),
      const std::set<TraversalOption> &allowedTraversalOptions = TraversalOption::getAllOptions(),
      const std::set<LoadEstimatorOption> &allowedLoadEstimatorOptions = LoadEstimatorOption::getAllOptions(),
      const std::set<DataLayoutOption> &allowedDataLayoutOptions = DataLayoutOption::getAllOptions(),
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(),
      const std::set<VectorizationPatternOption> &allowedVecPatternOptions =
          VectorizationPatternOption::getAllOptions(),
      size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      size_t predNumLHSamples = 1000, unsigned long seed = std::random_device()());

  TuningStrategyOption getOptionType() const override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  bool reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  bool optimizeSuggestions(std::vector<Configuration> &configQueue,
                           const EvidenceCollection &evidenceCollection) override;

  bool needsDomainSimilarityStatistics() const override;

  void rejectConfiguration(const Configuration &configuration, bool indefinitely) override;

  /**
   * Indicate if the search space contains any configurations.
   * @return
   */
  bool searchSpaceIsEmpty() const;

 private:
  /**
   * Generate n samples and sort them depending on
   * given acquistion function.
   * @param n numSamples
   * @param af acquisition function
   */
  std::vector<FeatureVector> sampleAcquisitions(size_t n, AcquisitionFunctionOption af);

  const InteractionTypeOption _interactionType;

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<VectorizationPatternOption> _vecPatternOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
  FeatureVectorEncoder _encoder;

  /**
   * Set of all configs that were marked as invalid.
   */
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;

  Random _rng;
  GaussianProcess _gaussianProcess;
  size_t _maxEvidence;
  AcquisitionFunctionOption _predAcqFunction;
  size_t _predNumLHSamples;
};

}  // namespace autopas
