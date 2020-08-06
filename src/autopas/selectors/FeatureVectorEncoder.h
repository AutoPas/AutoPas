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
   * Number of discrete values in a cluster encoded vector
   */
  static constexpr size_t clusterDiscreteSize = 3;
  /**
   * Index containing information about ContainerTraversalEstimator in continuous part of a cluster encoded vector
   */
  static constexpr size_t containerTraversalEstimatorClusterIndex = 0;
  /**
   * Index containing information about DataLayout in continuous part of a cluster encoded vector
   */
  static constexpr size_t dataLayoutClusterIndex = 1;
  /**
   * Index containing information about Newton3 in continuous part of a cluster encoded vector
   */
  static constexpr size_t newtonClusterIndex = 2;
  /**
   * Index containing information about CellSizeFactor in discrete part of a cluster encoded vector
   */
  static constexpr size_t cellSizeFactorClusterIndex = 0.;

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

    _dimRestrictions.clear();
    _dimRestrictions.reserve(clusterDiscreteSize);
    _dimRestrictions.push_back(_containerTraversalEstimatorOptions.size());
    _dimRestrictions.push_back(_dataLayoutOptions.size());
    _dimRestrictions.push_back(_newton3Options.size());
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
    auto [container, traversal, estimator] =
        _containerTraversalEstimatorOptions[vecDiscrete[containerTraversalEstimatorClusterIndex]];

    return FeatureVector(container, vecContinuous[cellSizeFactorClusterIndex], traversal, estimator,
                         _dataLayoutOptions[vecDiscrete[dataLayoutClusterIndex]],
                         _newton3Options[vecDiscrete[newtonClusterIndex]]);
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
    int containerTraversalEstimatorIndex = static_cast<int>(
        std::distance(_containerTraversalEstimatorOptions.begin(),
                      std::find(_containerTraversalEstimatorOptions.begin(), _containerTraversalEstimatorOptions.end(),
                                std::make_tuple(vec.container, vec.traversal, vec.loadEstimator))));
    int dataLayoutIndex = static_cast<int>(std::distance(
        _dataLayoutOptions.begin(), std::find(_dataLayoutOptions.begin(), _dataLayoutOptions.end(), vec.dataLayout)));
    int newton3Index = static_cast<int>(
        std::distance(_newton3Options.begin(), std::find(_newton3Options.begin(), _newton3Options.end(), vec.newton3)));

    Eigen::Vector3i vecDiscrete({containerTraversalEstimatorIndex, dataLayoutIndex, newton3Index});
    Eigen::VectorXd vecContinuous(FeatureVector::featureSpaceContinuousDims + 1);
    vecContinuous << vec.cellSizeFactor, iteration;
    return std::make_pair(vecDiscrete, vecContinuous);
  }

  /**
   * Inverse of convertToClusterWithIteration. Convert cluster representation with iteration back
   * to Feature vector while ignoring the iteration.
   * @param vec cluster encoded vector
   * @return decoded vector
   */
  FeatureVector convertFromClusterWithIteration(const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec) {
    const auto &[vecDiscrete, vecContinuous] = vec;
    auto [container, traversal, estimator] =
        _containerTraversalEstimatorOptions[vecDiscrete[containerTraversalEstimatorClusterIndex]];

    return FeatureVector(container, vecContinuous[cellSizeFactorClusterIndex], traversal, estimator,
                         _dataLayoutOptions[vecDiscrete[dataLayoutClusterIndex]],
                         _newton3Options[vecDiscrete[newtonClusterIndex]]);
  }

  /**
   * Get cluster-encoded neighbours of given target with fixed weight.
   * Neighbours are all configurations which differ in at most one configuration from target.
   * @param target
   * @return
   */
  std::vector<std::pair<Eigen::VectorXi, double>> clusterNeighboursManhattan1(const Eigen::VectorXi &target) {
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
  std::vector<std::pair<Eigen::VectorXi, double>> clusterNeighboursManhattan1Container(const Eigen::VectorXi &target) {
    std::vector<std::pair<Eigen::VectorXi, double>> result;
    // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
    // dimension)
    result.reserve(std::accumulate(_dimRestrictions.begin(), _dimRestrictions.end(), -_dimRestrictions.size()));

    auto targetContainer =
        std::get<0>(_containerTraversalEstimatorOptions[target[containerTraversalEstimatorClusterIndex]]);

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
          if (i == containerTraversalEstimatorClusterIndex) {
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
  std::vector<FeatureVector::ContainerTraversalEstimatorOption> _containerTraversalEstimatorOptions;
  std::vector<DataLayoutOption> _dataLayoutOptions;
  std::vector<Newton3Option> _newton3Options;

  /**
   * Restriction of a cluster-encoded vector in each dimension.
   */
  std::vector<int> _dimRestrictions;

  /**
   * Dimensions of a one-hot-encoded vector.
   */
  size_t _oneHotDims{0};
};

}  // namespace autopas
