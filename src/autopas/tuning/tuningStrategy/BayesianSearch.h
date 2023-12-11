/**
 * @file BayesianSearch.h
 * @author Jan Nguyen
 * @date 12.06.2019
 */

#pragma once

#include <set>
#include <unordered_set>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/GaussianModel/GaussianProcess.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyInterface.h"
#include "autopas/tuning/utils/FeatureVectorEncoder.h"
#include "autopas/utils/NumberInterval.h"
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
  constexpr static std::size_t maxAttempts = 10;
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
      const std::set<Newton3Option> &allowedNewton3Options = Newton3Option::getAllOptions(), std::size_t maxEvidence = 10,
      AcquisitionFunctionOption predAcqFunction = AcquisitionFunctionOption::upperConfidenceBound,
      std::size_t predNumLHSamples = 1000, unsigned long seed = std::random_device()());

  TuningStrategyOption getOptionType() override;

  void addEvidence(const Configuration &configuration, const Evidence &evidence) override;

  void reset(std::size_t iteration, std::size_t tuningPhase, std::vector<Configuration> &configQueue,
             const autopas::EvidenceCollection &evidenceCollection) override;

  void optimizeSuggestions(std::vector<Configuration> &configQueue, const EvidenceCollection &evidence) override;

  bool needsSmoothedHomogeneityAndMaxDensity() const override;

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
  std::vector<FeatureVector> sampleAcquisitions(std::size_t n, AcquisitionFunctionOption af);

  std::set<ContainerOption> _containerOptionsSet;
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;
  std::unique_ptr<NumberSet<double>> _cellSizeFactors;
  FeatureVectorEncoder _encoder;

  /**
   * Set of all configs that were marked as invalid.
   */
  std::unordered_set<FeatureVector, ConfigHash> _invalidConfigs;

  Random _rng;
  GaussianProcess _gaussianProcess;
  std::size_t _maxEvidence;
  AcquisitionFunctionOption _predAcqFunction;
  std::size_t _predNumLHSamples;
};

}  // namespace autopas
