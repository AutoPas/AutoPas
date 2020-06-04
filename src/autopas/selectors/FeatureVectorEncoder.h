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
 * Encoder to convert FeatureVector from and to Eigen:Vector
 */
class FeatureVectorEncoder {
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

    _oneHotDims = FeatureVector::featureSpaceContinuousDims + _containerTraversalEstimatorOptions.size() +
                  _dataLayoutOptions.size() + _newton3Options.size();
  }

  /**
   * Get the dimensions of a one-hot encoded vector.
   * @return
   */
  [[nodiscard]] size_t getOneHotDims() { return _oneHotDims; }

  /**
   * Encode FeatureVector to Eigen::VectorXd using one-hot-encoding.
   * @param vec vector to encode
   * @return one-hot-encoded vector
   */
  [[nodiscard]] Eigen::VectorXd oneHotEncode(const FeatureVector &vec) const {
    std::vector<double> data;
    data.reserve(_oneHotDims);

    data.push_back(vec.cellSizeFactor);
    for (const auto &[container, traversal, estimator] : _containerTraversalEstimatorOptions) {
      data.push_back(
          (container == vec.container and traversal == vec.traversal and estimator == vec.loadEstimator) ? 1. : 0.);
    }
    for (const auto &dataLayout : _dataLayoutOptions) {
      data.push_back((dataLayout == vec.dataLayout) ? 1. : 0.);
    }
    for (const auto &newton3 : _newton3Options) {
      data.push_back((newton3 == vec.newton3) ? 1. : 0.);
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
    double cellSizeFactor = vec[pos++];

    // get container, traversal, loadEstimator
    std::optional<FeatureVector::ContainerTraversalEstimatorOption> cteOption{};
    for (const auto &option : _containerTraversalEstimatorOptions) {
      if (vec[pos++] == 1.) {
        if (cteOption) {
          utils::ExceptionHandler::exception(
              "FeatureVectorEncoder.oneHotDecode: Vector encodes more than one ContainerTraversalEstimator. (More than "
              "one value for ContainerTraversalEstimator "
              "equals 1.)");
        }
        cteOption = option;
      }
    }
    if (not cteOption) {
      utils::ExceptionHandler::exception(
          "FeatureVectorEncoder.oneHotDecode: Vector encodes no ContainerTraversalEstimator. (All values for "
          "ContainerTraversalEstimator equal 0.)");
    }

    // get data layout
    std::optional<DataLayoutOption> dataLayout = {};
    for (const auto &option : _dataLayoutOptions) {
      if (vec[pos++] == 1.) {
        if (dataLayout) {
          utils::ExceptionHandler::exception(
              "FeatureVectorEncoder.oneHotDecode: Vector encodes more than one data layout. (More than one value for "
              "dataLayout equals 1.)");
        }
        dataLayout = option;
      }
    }
    if (not dataLayout) {
      utils::ExceptionHandler::exception(
          "FeatureVectorEncoder.oneHotDecode: Vector encodes no data layout. (All values for dataLayout equal 0.)");
    }

    // get newton3
    std::optional<Newton3Option> newton3 = {};
    for (const auto &option : _newton3Options) {
      if (vec[pos++] == 1.) {
        if (newton3) {
          utils::ExceptionHandler::exception(
              "FeatureVectorEncoder.oneHotDecode: Vector encodes more than one newton3. (More than one value for "
              "newton3 "
              "equals 1.)");
        }
        newton3 = option;
      }
    }
    if (not newton3) {
      utils::ExceptionHandler::exception(
          "FeatureVectorEncoder.oneHotDecode: Vector encodes no newton3. (All values for newton3 equal 0.)");
    }

    const auto &[container, traversal, estimator] = *cteOption;
    return FeatureVector(container, cellSizeFactor, traversal, estimator, *dataLayout, *newton3);
  }

  /**
   * Convert Feature vector to cluster representation for GaussianCluster.
   * Discrete values are encoded using their index in given std::vector.
   * @param vec vector to encode
   * @return cluster encoded vector
   */
  [[nodiscard]] std::pair<Eigen::VectorXi, Eigen::VectorXd> convertToCluster(const FeatureVector &vec) const {
    int containerTraversalEstimatorIndex = static_cast<int>(
        std::distance(_containerTraversalEstimatorOptions.begin(),
                      std::find(_containerTraversalEstimatorOptions.begin(), _containerTraversalEstimatorOptions.end(),
                                std::make_tuple(vec.container, vec.traversal, vec.loadEstimator))));
    int dataLayoutIndex = static_cast<int>(std::distance(
        _dataLayoutOptions.begin(), std::find(_dataLayoutOptions.begin(), _dataLayoutOptions.end(), vec.dataLayout)));
    int newton3Index = static_cast<int>(
        std::distance(_newton3Options.begin(), std::find(_newton3Options.begin(), _newton3Options.end(), vec.newton3)));

    Eigen::Vector3i vecDiscrete({containerTraversalEstimatorIndex, dataLayoutIndex, newton3Index});
    Eigen::VectorXd vecContinuous(FeatureVector::featureSpaceContinuousDims);
    vecContinuous << vec.cellSizeFactor;
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
    auto [container, traversal, estimator] = _containerTraversalEstimatorOptions[vecDiscrete[0]];

    return FeatureVector(container, vecContinuous[0], traversal, estimator, _dataLayoutOptions[vecDiscrete[1]],
                         _newton3Options[vecDiscrete[2]]);
  }

  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;

  /**
   * Dimensions of a one-hot-encoded vector
   */
  size_t _oneHotDims;
};

}  // namespace autopas
