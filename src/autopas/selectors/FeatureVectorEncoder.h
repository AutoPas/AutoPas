/**
 * @file FeatureVectorEncoder.h
 * @author Jan Nguyen
 * @date 04.06.20
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "FeatureVector.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Encoder to convert FeatureVector from and to Eigen::Vector.
 */
class FeatureVectorEncoder {
  /**
   * Number of tunable discrete dimensions.
   */
  static constexpr size_t tunableDiscreteDims = 3;
  /**
   * Index containing information about ContainerTraversalEstimator in discrete part of convertToTunable().
   */
  static constexpr size_t containerTraversalEstimatorIndex = 0;
  /**
   * Index containing information about DataLayout in discrete part of convertToTunable().
   */
  static constexpr size_t dataLayoutIndex = 1;
  /**
   * Index containing information about Newton3 in discrete part of convertToTunable(),
   */
  static constexpr size_t newtonIndex = 2;

  /**
   * Number of tunable continuous dimensions.
   */
  static constexpr size_t tunableContinuousDims = 1;
  /**
   * Index containing information about CellSizeFactor in continuous part of convertToTunable().
   */
  static constexpr size_t cellSizeFactorIndex = 0.;

  /**
   * Array which consists of an int for each tunable discrete dimension.
   */
  using DiscreteDimensionType = std::array<int, tunableDiscreteDims>;
  /**
   * Array which consists of a double for each tunable continuous dimension.
   */
  using ContinuousDimensionType = std::array<double, tunableContinuousDims>;
  /**
   * Pair of array for discrete and continuous dimensions.
   */
  using DiscreteContinuousPair = std::pair<DiscreteDimensionType, ContinuousDimensionType>;

 public:
  /**
   * Default Constructor
   */
  FeatureVectorEncoder() : _containerTraversalEstimatorOptions(), _dataLayoutOptions(), _newton3Options() {}

  /**
   * Contructor
   * @param containerTraversalEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   */
  FeatureVectorEncoder(
      const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options) {
    setAllowedOptions(containerTraversalEstimatorOptions, dataLayoutOptions, newton3Options);
  }

  /**
   * Set allowed options. All previously encoded vector can not be decoded anymore.
   * @param containerTraversalEstimatorOptions
   * @param dataLayoutOptions
   * @param newton3Options
   */
  void setAllowedOptions(
      const std::vector<FeatureVector::ContainerTraversalEstimatorOption> &containerTraversalEstimatorOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options) {
    _containerTraversalEstimatorOptions = containerTraversalEstimatorOptions;
    _dataLayoutOptions = dataLayoutOptions;
    _newton3Options = newton3Options;

    _oneHotDims = _containerTraversalEstimatorOptions.size() + _dataLayoutOptions.size() + _newton3Options.size() +
                  tunableContinuousDims;

    _dimRestrictions[containerTraversalEstimatorIndex] = _containerTraversalEstimatorOptions.size();
    _dimRestrictions[dataLayoutIndex] = _dataLayoutOptions.size();
    _dimRestrictions[newtonIndex] = _newton3Options.size();
  }

  /**
   * Get the dimensions of a one-hot encoded vector.
   * @return
   */
  [[nodiscard]] size_t getOneHotDims() const { return _oneHotDims; }

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
      for (int listIndex = 0; listIndex < _dimRestrictions[i]; ++listIndex) {
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
      for (int listIndex = 0; listIndex < _dimRestrictions[i]; ++listIndex) {
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
   * @param vec vector to encode
   * @return cluster encoded vector
   */
  [[nodiscard]] std::pair<Eigen::VectorXi, Eigen::VectorXd> convertToCluster(const FeatureVector &vec) const {
    auto [discreteValues, continuousValues] = convertToTunable(vec);
    Eigen::Map<Eigen::VectorXi> vecDiscrete(discreteValues.data(), tunableDiscreteDims);
    Eigen::Map<Eigen::VectorXd> vecContinuous(continuousValues.data(), tunableContinuousDims);
    return std::make_pair(vecDiscrete, vecContinuous);
  }

  /**
   * Inverse of convertToCluster. Convert cluster representation back
   * to Feature vector.
   * @param vec cluster encoded vector
   * @return decoded vector
   */
  FeatureVector convertFromCluster(const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec) {
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
   * Convert Feature vector to cluster representation for GaussianCluster.
   * Discrete values are encoded using their index in given std::vector.
   * Additionaly append current iteration to the continuous tuple.
   * @param vec vector to encode
   * @param iteration current iteration which may be scaled by some factor
   * @return cluster encoded vector
   */
  [[nodiscard]] std::pair<Eigen::VectorXi, Eigen::VectorXd> convertToClusterWithIteration(const FeatureVector &vec,
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
   * Inverse of convertToClusterWithIteration. Convert cluster representation with iteration back
   * to Feature vector while ignoring the iteration.
   * @param vec cluster encoded vector
   * @return decoded vector
   */
  [[nodiscard]] FeatureVector convertFromClusterWithIteration(const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec) {
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
    result.reserve(std::accumulate(_dimRestrictions.begin(), _dimRestrictions.end(), -_dimRestrictions.size()));

    // for each dimension
    for (int i = 0; i < target.size(); ++i) {
      // initial value
      auto init = target[i];

      // for each possible value of that dimension
      for (int x = 0; x < _dimRestrictions[i]; ++x) {
        // skip initial value
        if (x != init) {
          auto neighbour = target;
          neighbour[i] = x;
          result.push_back(std::make_pair(std::move(neighbour), 1.));
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
    result.reserve(std::accumulate(_dimRestrictions.begin(), _dimRestrictions.end(), -_dimRestrictions.size()));

    auto targetContainer = std::get<0>(_containerTraversalEstimatorOptions[target[containerTraversalEstimatorIndex]]);

    // for each dimension
    for (int i = 0; i < target.size(); ++i) {
      // initial value
      auto init = target[i];

      // for each possible value of that dimension
      for (int x = 0; x < _dimRestrictions[i]; ++x) {
        // skip initial value
        if (x != init) {
          auto neighbour = target;
          neighbour[i] = x;
          double weight = 1.;
          // check if container changed
          if (i == containerTraversalEstimatorIndex) {
            auto xContainer = std::get<0>(_containerTraversalEstimatorOptions[x]);
            if (targetContainer != xContainer) {
              // assign lower weight
              weight = 0.5;
            }
          }

          result.push_back(std::make_pair(std::move(neighbour), weight));
        }
      }
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
  [[nodiscard]] DiscreteContinuousPair convertToTunable(const FeatureVector &vec) const {
    DiscreteDimensionType discreteValues;
    discreteValues[containerTraversalEstimatorIndex] =
        getIndex(_containerTraversalEstimatorOptions, std::make_tuple(vec.container, vec.traversal, vec.loadEstimator));
    discreteValues[dataLayoutIndex] = getIndex(_dataLayoutOptions, vec.dataLayout);
    discreteValues[newtonIndex] = getIndex(_newton3Options, vec.newton3);

    ContinuousDimensionType continuousValues;
    continuousValues[cellSizeFactorIndex] = vec.cellSizeFactor;

    return std::make_pair(std::move(discreteValues), std::move(continuousValues));
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
        _containerTraversalEstimatorOptions[discreteValues[containerTraversalEstimatorIndex]];
    auto dataLayout = _dataLayoutOptions[discreteValues[dataLayoutIndex]];
    auto newton3 = _newton3Options[discreteValues[newtonIndex]];

    auto cellSizeFactor = continuousValues[cellSizeFactorIndex];

    return FeatureVector(container, cellSizeFactor, traversal, estimator, dataLayout, newton3);
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

  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;

  /**
   * Restriction of a cluster-encoded vector in each dimension.
   */
  std::array<int, tunableDiscreteDims> _dimRestrictions;

  /**
   * Dimensions of a one-hot-encoded vector.
   */
  size_t _oneHotDims{0};
};

}  // namespace autopas
