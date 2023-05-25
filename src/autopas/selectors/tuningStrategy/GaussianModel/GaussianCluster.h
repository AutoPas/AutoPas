/**
 * @file GaussianCluster.h
 * @author Jan Nguyen
 * @date 05.04.20
 */

#pragma once

#include <Eigen/Core>
#include <utility>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/selectors/tuningStrategy/GaussianModel/GaussianModelTypes.h"
#include "autopas/selectors/tuningStrategy/GaussianModel/GaussianProcess.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/logging/GaussianClusterLogger.h"

namespace autopas {

/**
 * Model to predicts the output of a blackbox function f(x) for given input x.
 * The model separates discrete and continuous dimensions of x. For each possible
 * discrete tuple a Gaussian process is assigned to estimate f(x) if the tuple is fixed.
 * Some sample input-output pairs (x,f(x)) should be provided as evidence.
 */
class GaussianCluster {
  // number of samples to find optimal hyperparameters
  static constexpr size_t hp_sample_size = 500;
  // number of hyperparameters
  static constexpr size_t hp_size = 25;

 public:
  /**
   * Different weight functions between clusters
   */
  enum WeightFunction {
    /**
     * geometric mean of probability densitiy over all evidence in neighbouring cluster if provided to the target
     * cluster
     */
    evidenceMatchingProbabilityGM,
    /**
     * geometric mean of scaled probability densitiy over all evidence in neighbouring cluster if provided to the target
     * cluster Probability densities are scaled such that the maximum is 1.
     */
    evidenceMatchingScaledProbabilityGM,
    /**
     * Wasserstein-2 distance of normal distributions of cluster given a continuous cluster
     */
    wasserstein2
  };

  /**
   * Constructor
   * @param dimRestriction restrict the i-th dimension to a integer between 0 and dimRestriction[i]-1
   * @param continuousDims additional unrestricted dimensions
   * @param weightFun function to calculate weight between clusters
   * @param sigma fixed noise
   * @param rngRef reference to random number generator
   * @param vectorToString function to convert vectors to a readable string
   * @param outputSuffix Suffix for all output files produced by this class.
   */
  GaussianCluster(const std::vector<int> &dimRestriction, size_t continuousDims, WeightFunction weightFun, double sigma,
                  Random &rngRef, const GaussianModelTypes::VectorToStringFun &vectorToString = defaultVecToString,
                  const std::string &outputSuffix = "");

  ~GaussianCluster();

  /**
   * Get the number of clusters in each dimension.
   * @return
   */
  [[nodiscard]] const std::vector<int> &getDimensions() const;

  /**
   * Change the number of cluster in all dimension.
   * This will discard all evidence.
   * @param newValue new number of clusters in each dimension
   */
  void setDimensions(const std::vector<int> &newValue);

  /**
   * Get the underlying GaussianProcess of a cluster.
   * This function should only be used for testing purposes.
   * @param index1D
   * @return
   */
  [[nodiscard]] const GaussianProcess &getCluster(size_t index1D) const;

  /**
   * Discard all evidence.
   */
  void clear();

  /**
   * Get the number of evidence provided.
   * @return
   */
  [[nodiscard]] size_t numEvidence() const;

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improve the quality of future predictions.
   * @param inputDiscrete x
   * @param inputContinuous y
   * @param output f((x,y))
   */
  void addEvidence(const GaussianModelTypes::VectorDiscrete &inputDiscrete,
                   const GaussianModelTypes::VectorContinuous &inputContinuous, double output);

  /**
   * Provide a input-output pair as evidence.
   * Each evidence improves the quality of future predictions.
   * @param input (x,y)
   * @param output f((x,y))
   */
  void addEvidence(const GaussianModelTypes::VectorPairDiscreteContinuous &input, double output);

  /**
   * Get the evidence with the highest output value
   * @return input of max
   */
  [[nodiscard]] GaussianModelTypes::VectorPairDiscreteContinuous getEvidenceMax() const;

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples
   * and calculate their corresponding acquisition.
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param continuousSamples continuous tuples
   * @return all discrete-continuous tuples paired with their corresponding acquisition
   */
  [[nodiscard]] std::vector<GaussianModelTypes::VectorAcquisition> sampleAcquisition(
      AcquisitionFunctionOption af, const GaussianModelTypes::NeighbourFunction &neighbourFun,
      const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const;

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples
   * and calculate their weight to each other. Output graph in AutoPasLog Trace.
   * @param neighbourFun function which generates neighbours of given discrete tuple
   * @param continuousSamples continuous tuples
   */
  void logDebugGraph(const GaussianModelTypes::NeighbourFunction &neighbourFun,
                     const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const;
  /**
   * Change the used function to convert from vector to string.
   * @param fun new converter
   */
  void setVectorToStringFun(const GaussianModelTypes::VectorToStringFun &fun);

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples and
   * returns the vector with maximum acquisition.
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours with prior weight of given discrete tuple
   * @param continuousSamples continuous tuples
   * @return pair of vector and corresponding maximum acquistion
   */
  [[nodiscard]] GaussianModelTypes::VectorAcquisition sampleAcquisitionMax(
      AcquisitionFunctionOption af, const GaussianModelTypes::NeighbourFunction &neighbourFun,
      const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const;

  /**
   * Generate all possible combinations of discrete tuples and continuous tuples in samples and
   * order them by their acquisition.
   * @param af function to optimize
   * @param neighbourFun function which generates neighbours with prior weight of given discrete tuple
   * @param continuousSamples continuous tuples
   * @return all discrete-continuous tuples ordered by their corresponding acquisition
   */
  [[nodiscard]] std::vector<GaussianModelTypes::VectorPairDiscreteContinuous> sampleOrderedByAcquisition(
      AcquisitionFunctionOption af, const GaussianModelTypes::NeighbourFunction &neighbourFun,
      const std::vector<GaussianModelTypes::VectorContinuous> &continuousSamples) const;

  /**
   * Default function used to convert vectors to readable strings.
   * @note begin() and end() currently not available for Eigen::Vector, so AutoPas ArrayUtils cannot be used.
   * @param vec
   * @return string with format (a,b,...,n) beginning with discrete values.
   */
  static std::string defaultVecToString(const GaussianModelTypes::VectorPairDiscreteContinuous &vec);

 private:
  /**
   * Create a GaussianProcess for each cluster and precalculate DiscreteVector for each cluster index.
   */
  void initClusters();

  /**
   * Get the cluster index of a discrete tuple.
   * @param x the discrete tuple
   * @return
   */
  [[nodiscard]] size_t getIndex(const GaussianModelTypes::VectorDiscrete &x) const;

  /**
   * Increment the given discrete tuple x, resulting
   * in the tuple following x in the cluster list.
   * The increment works in a date-like fashion:
   * One increment adds one to the first vector entry.
   * If this entry overflows also the next entry is incremented and so on.
   * @param x discrete tuple
   */
  void discreteIncrement(GaussianModelTypes::VectorDiscrete &x) const;

  /**
   * Calculate mean, variance and stddev for all clusters for given continous tuple.
   * @param continuousTuple
   * @return Vectors means, vars, stddevs containing corresponding values for each cluster
   */
  [[nodiscard]] std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> precalculateDistributions(
      const GaussianModelTypes::VectorContinuous &continuousTuple) const;

  /**
   * Calculate acquisition for all clusters.
   * @param af acquisition function
   * @param means mean for each cluster
   * @param vars variance for each cluster
   * @return
   */
  [[nodiscard]] std::vector<double> precalculateAcquisitions(AcquisitionFunctionOption af,
                                                             const std::vector<double> &means,
                                                             const std::vector<double> &vars) const;

  /**
   * Initalize the neighbour-weight list for each cluster.
   * @param neighbourFun function which generates neighbours with prior weight of given discrete tuple
   * @return
   */
  [[nodiscard]] GaussianModelTypes::NeighboursWeights initNeighbourWeights(
      const GaussianModelTypes::NeighbourFunction &neighbourFun) const;

  /**
   * Update the neighbour-weight list for each cluster. This function is called for each continuous tuple.
   * Given mean, var and stddev should be evaluated for the current continuous tuple.
   * @param neighbourWeights list to update
   * @param means mean for each cluster
   * @param vars variance for each cluster
   * @param stddevs standard deviation for each cluster
   */
  void updateNeighbourWeights(GaussianModelTypes::NeighboursWeights &neighbourWeights, const std::vector<double> &means,
                              const std::vector<double> &vars, const std::vector<double> &stddevs) const;

  /**
   * Number of clusters per discrete dimension.
   */
  std::vector<int> _dimRestriction;
  /**
   * Number of additional unrestricted continuous dimensions.
   */
  const size_t _continuousDims;

  /**
   * Gaussian process for each discrete tuple.
   */
  std::vector<GaussianProcess> _clusters;

  /**
   * Vector to easiliy map from index to DiscreteVector
   */
  std::vector<GaussianModelTypes::VectorDiscrete> _discreteVectorMap;

  /**
   * Function used to calculate the weight between two clusters.
   */
  WeightFunction _weightFun;

  /**
   * Current smallest evidence output.
   */
  double _evidenceMinValue;
  /**
   * Current greatest evidence output.
   */
  double _evidenceMaxValue;
  /**
   * Current greatest evidence input
   */
  GaussianModelTypes::VectorPairDiscreteContinuous _evidenceMaxVector;
  /**
   * Current number of evidence.
   */
  size_t _numEvidence;

  /**
   * Fixed noise assumed.
   */
  const double _sigma;

  Random &_rng;

  /**
   * Logger for graphs. Can be used for visualization of cluster map.
   */
  std::unique_ptr<GaussianClusterLogger> _logger;
};
}  // namespace autopas
