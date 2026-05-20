/**
 * @file FeatureVectorEncoder.h
 * @author Jan Nguyen
 * @date 04.06.20
 */

#pragma once

#include <Eigen/Core>
#include <numeric>
#include <optional>
#include <vector>

#include "FeatureVector.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Encoder to convert FeatureVector from and to Eigen::Vector.
 */
class FeatureVectorEncoder {
  /**
   * One-hot encoded indices of convertToTunable(). These are categorical components of the configuration which
   * should not be treated as ordinal values (e.g., mapping Linked Cell C01 to 0.0, C08 to 1.0, and Verlet Lists
   * to 2.0 would be semantically incorrect). Instead, one-hot encoding is used: for each categorical option,
   * there is one binary dimension per possible value—exactly one dimension is 1 if the FeatureVector represents
   * that choice, and all others are 0.
   */
  enum class OneHotEncodedIndices { containerTraversalEstimator, dataLayout, newton3, ompKind, TOTALNUMBER };

  /**
   * Indices of the continuous part of convertToTunable().
   * Does not have to represent a continuous variable, but rather something where it makes sense to map to a continuous
   * variable rather than using one-hot encoding.
   */
  enum class ContinuousIndices { cellSizeFactor, ompChunkSize, TOTALNUMBER };

 public:
  /**
   * Number of tunable one-hot dimensions.
   */
  static constexpr size_t numTunableOneHotIndices{static_cast<size_t>(OneHotEncodedIndices::TOTALNUMBER)};

  /**
   * Number of tunable continuous dimensions.
   */
  static constexpr size_t numTunableContinuousDims{static_cast<size_t>(ContinuousIndices::TOTALNUMBER)};

 private:
  /**
   * Generalized form of all one-hot vector representations. Every representation of one-hot values should be able to
   * be brought to and from this form. Note: this is not the one-hot form, just simply a form where integers are used
   * instead of e.g. enums, for each index.
   */
  using OneHotIndicesType = std::array<int, numTunableOneHotIndices>;
  /**
   * Generalized form of all continuous vector representations. Any representation of continuous values should be able
   * to be brought to and from this form.
   */
  using ContinuousIndicesType = std::array<double, numTunableContinuousDims>;

 public:
  /**
   * Default Constructor
   */
  FeatureVectorEncoder();

  /**
   * Constructor
   * @param containerTraversalEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   * @param cellSizeFactors
   * @param ompKindOptions
   * @param ompChunkSizes
   * @param interactionType
   */
  FeatureVectorEncoder(
      const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
      const NumberSet<double> &cellSizeFactors, const std::vector<OpenMPKindOption> &ompKindOptions,
      const NumberSet<size_t> &ompChunkSizes, const InteractionTypeOption &interactionType);

  ~FeatureVectorEncoder();

  /**
   * Set allowed options. All previously encoded vector can not be decoded anymore.
   * @param containerTraversalEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   * @param cellSizeFactors
   * @param ompKindOptions
   * @param ompChunkSizes
   */
  void setAllowedOptions(
      const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
      const NumberSet<double> &cellSizeFactors, const std::vector<OpenMPKindOption> &ompKindOptions,
      const NumberSet<size_t> &ompChunkSizes);

  /**
   * Get the dimensions of a one-hot encoded vector.
   * @return
   */
  [[nodiscard]] size_t getOneHotDims() const;

  /**
   * Get the number of allowed options of each discrete dimension.
   * @return
   */
  [[nodiscard]] const std::array<int, numTunableOneHotIndices> &getDiscreteRestrictions() const;

  /**
   * Encode FeatureVector to Eigen::VectorXd using one-hot-encoding.
   * @param vec vector to encode
   * @return one-hot-encoded vector
   */
  [[nodiscard]] Eigen::VectorXd oneHotEncode(const FeatureVector &vec) const;

  /**
   * Decode one-hot-encoded VectorXd to FeatureVector.
   * @param vec one-hot-encoded vector
   * @return decoded FeatureVector
   */
  [[nodiscard]] FeatureVector oneHotDecode(const Eigen::VectorXd &vec);

  /**
   * Convert Feature vector to cluster representation for GaussianCluster.
   * Discrete values are encoded using their index in given std::vector.
   * Additionally, append current iteration to the continuous tuple.
   * @param vec vector to encode
   * @param iteration current iteration which may be scaled by some factor
   * @return cluster encoded vector
   */
  [[nodiscard]] std::pair<Eigen::VectorXi, Eigen::VectorXd> convertToCluster(const FeatureVector &vec,
                                                                             double iteration) const;

  /**
   * Inverse of convertToCluster. Convert cluster representation back
   * to Feature vector while ignoring the iteration.
   * @param vec cluster encoded vector
   * @return decoded vector
   */
  [[nodiscard]] FeatureVector convertFromCluster(const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec);

  /**
   * Get cluster-encoded neighbors of given target with fixed weight.
   * Neighbors are all configurations which differ in at most one configuration from target.
   * @param target
   * @return
   */
  [[nodiscard]] std::vector<std::pair<Eigen::VectorXi, double>> clusterNeighborsManhattan1(
      const Eigen::VectorXi &target);

  /**
   * Get cluster-encoded neighbors of given target. Neighbors are all configurations which differ
   * in at most one configuration from target. The weight is lowered if container is changed.
   * @param target
   * @return
   */
  [[nodiscard]] std::vector<std::pair<Eigen::VectorXi, double>> clusterNeighborsManhattan1Container(
      const Eigen::VectorXi &target);

  /**
   * Create n latin-hypercube-samples from given featureSpace.
   * @param n number of samples
   * @param rng
   * @return vector of sample featureVectors
   */
  std::vector<FeatureVector> lhsSampleFeatures(size_t n, Random &rng) const;

  /**
   * Create n latin-hypercube-samples from the continuous featureSpace and append a value representing the
   * current iteration to each sample.
   * @param n number of samples
   * @param rng
   * @param iteration Current iteration which may be scaled by some factor.
   * @return vector of continuous feature samples
   */
  std::vector<Eigen::VectorXd> lhsSampleFeatureCluster(size_t n, Random &rng, double iteration) const;

 private:
  /**
   * Create a pair of array from a FeatureVector.
   * The first array contains all discrete options converted to a int which
   * indicates the position of the option in the respectively provided allowedOptions list.
   * The second array contains all continuous options which are converted to a double.
   * @param vec
   * @return
   */
  [[nodiscard]] std::pair<OneHotIndicesType, ContinuousIndicesType> convertToTunable(
      const FeatureVector &vec) const;

  /**
   * Inverse of convertToTunable. Generate FeatureVector from two arrays
   * which contain information of all discrete options and continuous options respectively.
   * @param discreteValues
   * @param continuousValues
   * @return
   */
  [[nodiscard]] FeatureVector convertFromTunable(const OneHotIndicesType &discreteValues,
                                                 const ContinuousIndicesType &continuousValues) const;

  /**
   * Get position of a value in provided list.
   * @param list
   * @param value
   * @return index i such that list[i] == value
   */
  template <typename T>
  static size_t getIndex(const std::vector<T> &list, const T &value) {
    for (size_t i = 0; i < list.size(); ++i) {
      if (list[i] == value) return i;
    }

    utils::ExceptionHandler::exception("FeatureVectorEncoder.getIndex: Value not allowed!");
    return list.size();
  }

  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions{};
  std::vector<DataLayoutOption> _dataLayoutOptions{};
  std::vector<Newton3Option> _newton3Options{};
  std::vector<OpenMPKindOption> _openMPKindOptions{};
  InteractionTypeOption _interactionType;
  /**
   * Number of allowed options of each discrete dimension.
   */
  std::array<int, numTunableOneHotIndices> _discreteRestrictions{};

  /**
   * Allowed values of each continuous dimension.
   */
  std::array<std::unique_ptr<NumberSet<double>>, numTunableContinuousDims> _continuousRestrictions{};

  /**
   * Dimensions of a one-hot-encoded vector.
   */
  size_t _oneHotDims{0};
};

}  // namespace autopas
