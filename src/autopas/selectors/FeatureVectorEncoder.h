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
   * Indices of the discrete part of convertToTunable().
   */
  enum class DiscreteIndices { containerTraversalEstimator, dataLayout, newton3, verletRebuildFrequency, TOTALNUMBER};

  /**
   * Indices of the continuous part of convertToTunable().
   */
  enum class ContinuousIndices { cellSizeFactor, TOTALNUMBER};

 public:
  /**
   * Number of tunable discrete dimensions.
   */
  static constexpr size_t tunableDiscreteDims{static_cast<size_t>(DiscreteIndices::TOTALNUMBER)};

  /**
   * Number of tunable continuous dimensions.
   */
  static constexpr size_t tunableContinuousDims{static_cast<size_t>(ContinuousIndices::TOTALNUMBER)};

 private:
  /**
   * Generalized form of all discrete vector representations. Every representation of discrete values should be able to
   * be brought to and from this form.
   */
  using DiscreteDimensionType = std::array<int, tunableDiscreteDims>;
  /**
   * Generalized form of all continuous vector representations. Any representation of continuous values should be able
   * to be brought to and from this form.
   */
  using ContinuousDimensionType = std::array<double, tunableContinuousDims>;

 public:
  /**
   * Default Constructor
   */
  FeatureVectorEncoder() {}

  /**
   * Contructor
   * @param containerTraversalEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   * @param cellSizeFactors
   * @param verletRebuildFrequencies
   */
  FeatureVectorEncoder(
      const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
      const NumberSet<double> &cellSizeFactors, const NumberSet<int> &verletRebuildFrequencies) {
    setAllowedOptions(containerTraversalEstimatorOptions, dataLayoutOptions, newton3Options, cellSizeFactors, verletRebuildFrequencies);
  }

  /**
   * Set allowed options. All previously encoded vector can not be decoded anymore.
   * @param containerTraversalEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   * @param cellSizeFactors
   * @param verletRebuildFrequencies
   */
  void setAllowedOptions(
      const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options,
      const NumberSet<double> &cellSizeFactors, const NumberSet<int> &verletRebuildFrequencies) {
    _containerTraversalEstimatorOptions = containerTraversalEstimatorOptions;
    _dataLayoutOptions = dataLayoutOptions;
    _newton3Options = newton3Options;
    auto verletRebuildFrequenciestmp = verletRebuildFrequencies.getAll();
    _verletRebuildFrequencies = std::vector<int>(verletRebuildFrequenciestmp.begin(), verletRebuildFrequenciestmp.end());

    _oneHotDims = _containerTraversalEstimatorOptions.size() + _dataLayoutOptions.size() + _newton3Options.size() +
                  _verletRebuildFrequencies.size() + tunableContinuousDims;

    _discreteRestrictions[static_cast<size_t>(DiscreteIndices::containerTraversalEstimator)] =
        _containerTraversalEstimatorOptions.size();
    _discreteRestrictions[static_cast<size_t>(DiscreteIndices::dataLayout)] = _dataLayoutOptions.size();
    _discreteRestrictions[static_cast<size_t>(DiscreteIndices::newton3)] = _newton3Options.size();

    _continuousRestrictions[static_cast<size_t>(ContinuousIndices::cellSizeFactor)] = cellSizeFactors.clone();
    _discreteRestrictions[static_cast<size_t>(DiscreteIndices::verletRebuildFrequency)] = verletRebuildFrequencies.size();
  }

  /**
   * Get the dimensions of a one-hot encoded vector.
   * @return
   */
  [[nodiscard]] size_t getOneHotDims() const { return _oneHotDims; }

  /**
   * Get the number of allowed options of each discrete dimension.
   * @return
   */
  [[nodiscard]] const std::array<int, tunableDiscreteDims> getDiscreteRestrictions() const {
    return _discreteRestrictions;
  }

  /**
   * Encode FeatureVector to Eigen::VectorXd using one-hot-encoding.
   * @param vec vector to encode
   * @return one-hot-encoded vector
   */
  [[nodiscard]] Eigen::VectorXd oneHotEncode(const FeatureVector &vec) const {
    std::vector<double> data;
    data.reserve(_oneHotDims);

    auto [discreteValues, continuousValues] = convertToTunable(vec);

    // discrete values are encoded using one-hot-encoding
    for (size_t i = 0; i < tunableDiscreteDims; ++i) {
      for (int listIndex = 0; listIndex < _discreteRestrictions[i]; ++listIndex) {
        data.push_back(listIndex == discreteValues[i] ? 1. : 0.);
      }
    }

    // continuous values are simply copied
    for (size_t i = 0; i < tunableContinuousDims; ++i) {
      data.push_back(continuousValues[i]);
    }

    return Eigen::Map<Eigen::VectorXd>(data.data(), _oneHotDims);
  }

  /**
   * Decode one-hot-encoded VectorXd to FeatureVector.
   * @param vec one-hot-encoded vector
   * @return decoded FeatureVector
   */
  [[nodiscard]] FeatureVector oneHotDecode(const Eigen::VectorXd &vec) {
    if (static_cast<size_t>(vec.size()) != _oneHotDims) {
      utils::ExceptionHandler::exception("FeatureVectorEncoder.oneHotDecode: Expected size {}, got {}", _oneHotDims,
                                         vec.size());
    }

    size_t pos = 0;

    DiscreteDimensionType discreteValues;

    // extract each one-hot-encoded discrete option
    for (size_t i = 0; i < tunableDiscreteDims; ++i) {
      // for each dimension get index whose value equals 1.
      std::optional<int> value = {};
      for (int listIndex = 0; listIndex < _discreteRestrictions[i]; ++listIndex) {
        if (vec[pos++] == 1.) {
          if (value) {
            utils::ExceptionHandler::exception(
                "FeatureVectorEncoder.oneHotDecode: Vector encodes more than one option for tunable dimension {}. "
                "(More than one value equals 1.)",
                i);
          }
          value = listIndex;
        }
      }
      if (not value) {
        utils::ExceptionHandler::exception(
            "FeatureVectorEncoder.oneHotDecode: Vector encodes no option for tunable dimension {}. (All values equal "
            "0.)",
            i);
      }
      discreteValues[i] = *value;
    }

    ContinuousDimensionType continuousValues;
    for (size_t i = 0; i < tunableContinuousDims; ++i) {
      continuousValues[i] = vec[pos++];
    }

    return convertFromTunable(discreteValues, continuousValues);
  }

  /**
   * Convert Feature vector to cluster representation for GaussianCluster.
   * Discrete values are encoded using their index in given std::vector.
   * Additionally, append current iteration to the continuous tuple.
   * @param vec vector to encode
   * @param iteration current iteration which may be scaled by some factor
   * @return cluster encoded vector
   */
  [[nodiscard]] std::pair<Eigen::VectorXi, Eigen::VectorXd> convertToCluster(const FeatureVector &vec,
                                                                             double iteration) const {
    auto [discreteValues, continuousValues] = convertToTunable(vec);
    Eigen::Map<Eigen::VectorXi> vecDiscrete(discreteValues.data(), tunableDiscreteDims);

    std::vector<double> continuousData;
    continuousData.reserve(tunableContinuousDims + 1);
    for (size_t i = 0; i < tunableContinuousDims; ++i) {
      continuousData.push_back(continuousValues[i]);
    }
    continuousData.push_back(iteration);
    Eigen::Map<Eigen::VectorXd> vecContinuous(continuousData.data(), tunableContinuousDims + 1);
    return std::make_pair(vecDiscrete, vecContinuous);
  }

  /**
   * Inverse of convertToCluster. Convert cluster representation back
   * to Feature vector while ignoring the iteration.
   * @param vec cluster encoded vector
   * @return decoded vector
   */
  [[nodiscard]] FeatureVector convertFromCluster(const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec) {
    const auto &[vecDiscrete, vecContinuous] = vec;

    DiscreteDimensionType discreteValues;
    for (size_t i = 0; i < tunableDiscreteDims; ++i) {
      discreteValues[i] = vecDiscrete[i];
    }

    ContinuousDimensionType continuousValues;
    for (size_t i = 0; i < tunableContinuousDims; ++i) {
      continuousValues[i] = vecContinuous[i];
    }

    return convertFromTunable(discreteValues, continuousValues);
  }

  /**
   * Get cluster-encoded neighbours of given target with fixed weight.
   * Neighbours are all configurations which differ in at most one configuration from target.
   * @param target
   * @return
   */
  [[nodiscard]] std::vector<std::pair<Eigen::VectorXi, double>> clusterNeighboursManhattan1(
      const Eigen::VectorXi &target) {
    std::vector<std::pair<Eigen::VectorXi, double>> result;
    // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
    // dimension)
    result.reserve(
        std::accumulate(_discreteRestrictions.begin(), _discreteRestrictions.end(), -_discreteRestrictions.size()));

    // for each dimension
    for (int i = 0; i < target.size(); ++i) {
      // initial value
      auto init = target[i];

      // for each possible value of that dimension
      for (int x = 0; x < _discreteRestrictions[i]; ++x) {
        // skip initial value
        if (x != init) {
          auto neighbour = target;
          neighbour[i] = x;
          result.emplace_back(std::move(neighbour), 1.);
        }
      }
    }
    return result;
  }

  /**
   * Get cluster-encoded neighbours of given target. Neighbours are all configurations which differ
   * in at most one configuration from target. The weight is lowered if container is changed.
   * @param target
   * @return
   */
  [[nodiscard]] std::vector<std::pair<Eigen::VectorXi, double>> clusterNeighboursManhattan1Container(
      const Eigen::VectorXi &target) {
    std::vector<std::pair<Eigen::VectorXi, double>> result;
    // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
    // dimension)
    result.reserve(
        std::accumulate(_discreteRestrictions.begin(), _discreteRestrictions.end(), -_discreteRestrictions.size()));

    auto targetContainer = std::get<0>(
        _containerTraversalEstimatorOptions[target[static_cast<int>(DiscreteIndices::containerTraversalEstimator)]]);

    // for each dimension
    for (int i = 0; i < target.size(); ++i) {
      // initial value
      auto init = target[i];

      // for each possible value of that dimension
      for (int x = 0; x < _discreteRestrictions[i]; ++x) {
        // skip initial value
        if (x != init) {
          auto neighbour = target;
          neighbour[i] = x;
          double weight = 1.;
          // check if container changed
          if (i == static_cast<int>(DiscreteIndices::containerTraversalEstimator)) {
            auto xContainer = std::get<0>(_containerTraversalEstimatorOptions[x]);
            if (targetContainer != xContainer) {
              // assign lower weight
              weight = 0.5;
            }
          }

          result.emplace_back(std::move(neighbour), weight);
        }
      }
    }
    return result;
  }

  /**
   * Create n latin-hypercube-samples from given featureSpace.
   * @param n number of samples
   * @param rng
   * @return vector of sample featureVectors
   */
  std::vector<FeatureVector> lhsSampleFeatures(size_t n, Random &rng) const {
    // create n samples for each discrete dimension.
    std::array<std::vector<size_t>, tunableDiscreteDims> lhsDiscreteSamples;
    for (size_t i = 0; i < tunableDiscreteDims; ++i) {
      lhsDiscreteSamples[i] = rng.uniformSample(0, _discreteRestrictions[i] - 1, n);
    }

    // create n samples for each continuous dimension.
    std::array<std::vector<double>, tunableContinuousDims> lhsContinuousSamples;
    for (size_t i = 0; i < tunableContinuousDims; ++i) {
      lhsContinuousSamples[i] = _continuousRestrictions[i]->uniformSample(n, rng);
    }

    // create FeatureVectors from raw samples
    std::vector<FeatureVector> result;
    for (size_t i = 0; i < n; ++i) {
      DiscreteDimensionType discreteValues;
      for (size_t d = 0; d < tunableDiscreteDims; ++d) {
        discreteValues[d] = lhsDiscreteSamples[d][i];
      }

      ContinuousDimensionType continuousValues;
      for (size_t c = 0; c < tunableContinuousDims; ++c) {
        continuousValues[c] = lhsContinuousSamples[c][i];
      }

      result.emplace_back(convertFromTunable(discreteValues, continuousValues));
    }

    return result;
  }

  /**
   * Create n latin-hypercube-samples from the continuous featureSpace and append a value representing the
   * current iteration to each sample.
   * @param n number of samples
   * @param rng
   * @param iteration Current iteration which may be scaled by some factor.
   * @return vector of continuous feature samples
   */
  std::vector<Eigen::VectorXd> lhsSampleFeatureCluster(size_t n, Random &rng, double iteration) const {
    // create n samples for each continuous dimension.
    std::array<std::vector<double>, tunableContinuousDims> lhsContinuousSamples;
    for (size_t i = 0; i < tunableContinuousDims; ++i) {
      lhsContinuousSamples[i] = _continuousRestrictions[i]->uniformSample(n, rng);
    }

    // create FeatureVectors from raw samples
    std::vector<Eigen::VectorXd> result;
    for (size_t i = 0; i < n; ++i) {
      std::array<double, tunableContinuousDims + 1> data;
      for (size_t c = 0; c < tunableContinuousDims; ++c) {
        data[c] = lhsContinuousSamples[c][i];
      }
      data[tunableContinuousDims] = iteration;
      result.emplace_back(Eigen::Map<Eigen::VectorXd>(data.data(), data.size()));
    }

    return result;
  }

 private:
  /**
   * Create a pair of array from a FeatureVector.
   * The first array contains all discrete options converted to a int which
   * indicates the position of the option in the respectively provided allowedOptions list.
   * The second array contains all continuous options which are converted to a double.
   * @param vec
   * @return
   */
  [[nodiscard]] std::pair<DiscreteDimensionType, ContinuousDimensionType> convertToTunable(
      const FeatureVector &vec) const {
    DiscreteDimensionType discreteValues;
    discreteValues[static_cast<size_t>(DiscreteIndices::containerTraversalEstimator)] =
        getIndex(_containerTraversalEstimatorOptions, std::make_tuple(vec.container, vec.traversal, vec.loadEstimator));
    discreteValues[static_cast<size_t>(DiscreteIndices::dataLayout)] = getIndex(_dataLayoutOptions, vec.dataLayout);
    discreteValues[static_cast<size_t>(DiscreteIndices::newton3)] = getIndex(_newton3Options, vec.newton3);

    ContinuousDimensionType continuousValues;
    continuousValues[static_cast<size_t>(ContinuousIndices::cellSizeFactor)] = vec.cellSizeFactor;
    discreteValues[static_cast<size_t>(DiscreteIndices::verletRebuildFrequency)] = getIndex(_verletRebuildFrequencies, vec.verletRebuildFrequency);

    return std::make_pair(discreteValues, continuousValues);
  }

  /**
   * Inverse of convertToTunable. Generate FeatureVector from two arrays
   * which contain information of all discrete options and continuous options respectively.
   * @param discreteValues
   * @param continuousValues
   * @return
   */
  [[nodiscard]] FeatureVector convertFromTunable(const DiscreteDimensionType &discreteValues,
                                                 const ContinuousDimensionType &continuousValues) const {
    const auto &[container, traversal, estimator] =
        _containerTraversalEstimatorOptions[discreteValues[static_cast<size_t>(
            DiscreteIndices::containerTraversalEstimator)]];
    auto dataLayout = _dataLayoutOptions[discreteValues[static_cast<size_t>(DiscreteIndices::dataLayout)]];
    auto newton3 = _newton3Options[discreteValues[static_cast<size_t>(DiscreteIndices::newton3)]];

    auto cellSizeFactor = continuousValues[static_cast<size_t>(ContinuousIndices::cellSizeFactor)];
    auto verletRebuildFrequency = _verletRebuildFrequencies[discreteValues[static_cast<size_t>(DiscreteIndices::verletRebuildFrequency)]];

    return FeatureVector(container, cellSizeFactor, traversal, estimator, dataLayout, newton3, verletRebuildFrequency);
  }

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
  std::vector<int> _verletRebuildFrequencies{};

  /**
   * Number of allowed options of each discrete dimension.
   */
  std::array<int, tunableDiscreteDims> _discreteRestrictions{};

  /**
   * Allowed values of each continuous dimension.
   */
  std::array<std::unique_ptr<NumberSet<double>>, tunableContinuousDims> _continuousRestrictions{};

  /**
   * Dimensions of a one-hot-encoded vector.
   */
  size_t _oneHotDims{0};
};

}  // namespace autopas
