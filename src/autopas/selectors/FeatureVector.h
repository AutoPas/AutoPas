/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "autopas/selectors/Configuration.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"

namespace autopas {

/**
 * FeatureVector representation of a Configuration
 */
class FeatureVector : public Configuration {
 public:
  /**
   * Number of tune-able dimensions.
   */
  static constexpr size_t featureSpaceDims = 5;

  /**
   * Dimensions of a one-hot-encoded vector
   * = 1 (cellSizeFactor) + traversals + loadEstimators + dataLayouts + newton3
   */
  inline static size_t oneHotDims = 1 + TraversalOption::getOptionNames().size() +
                                    LoadEstimatorOption::getOptionNames().size() +
                                    DataLayoutOption::getOptionNames().size() + Newton3Option::getOptionNames().size();

  /**
   * Default constructor. Results in invalid vector.
   */
  FeatureVector() : Configuration() {}

  /**
   * Constructor
   * @param container
   * @param traversal
   * @param loadEstimator
   * @param dataLayout
   * @param newton3
   * @param cellSizeFactor
   */
  FeatureVector(ContainerOption container, double cellSizeFactor, TraversalOption traversal,
                LoadEstimatorOption loadEstimator, DataLayoutOption dataLayout, Newton3Option newton3)
      : Configuration(container, cellSizeFactor, traversal, loadEstimator, dataLayout, newton3) {}

  /**
   * Construct from Configuration.
   * @param conf
   */
  FeatureVector(Configuration conf) : Configuration(conf) {}

  /**
   * Distance between two FeatureVectors.
   * Since there is no real ordering all discrete options are assumed to have a distance
   * of one to each other.
   * This function ignores the container dimension since it is encoded in the traversal.
   * @param other
   * @return
   */
  Eigen::VectorXd operator-(const FeatureVector &other) const {
    Eigen::VectorXd result(featureSpaceDims);
    result << cellSizeFactor - other.cellSizeFactor, traversal == other.traversal ? 0. : 1.,
        loadEstimator == other.loadEstimator ? 0. : 1., dataLayout == other.dataLayout ? 0. : 1.,
        newton3 == other.newton3 ? 0. : 1.;

    return result;
  }

  /**
   * Cast to Eigen::VectorXd ignoring ContainerOption.
   * @return
   */
  operator Eigen::VectorXd() const {
    Eigen::VectorXd result(featureSpaceDims);
    result << cellSizeFactor, static_cast<double>(traversal), static_cast<double>(loadEstimator),
        static_cast<double>(dataLayout), static_cast<double>(newton3);

    return result;
  }

  /**
   * Encode to Eigen::VectorXd ignoring ContainerOption using one-hot-encoding.
   * @return one-hot-encoded vector
   */
  Eigen::VectorXd oneHotEncode() const {
    std::vector<double> data;
    data.reserve(oneHotDims);

    data.push_back(cellSizeFactor);
    for (auto &[option, _] : TraversalOption::getOptionNames()) {
      data.push_back((option == traversal) ? 1. : 0.);
    }
    for (auto &[option, _] : LoadEstimatorOption::getOptionNames()) {
      data.push_back((option == loadEstimator) ? 1. : 0.);
    }
    for (auto &[option, _] : DataLayoutOption::getOptionNames()) {
      data.push_back((option == dataLayout) ? 1. : 0.);
    }
    for (auto &[option, _] : Newton3Option::getOptionNames()) {
      data.push_back((option == newton3) ? 1. : 0.);
    }

    return Eigen::Map<Eigen::VectorXd>(data.data(), oneHotDims);
  }

  /**
   * Decode one-hot-encoded VectorXd to FeatureVector.
   * Encoding ignores ContainerOption and valid options are unknown.
   * So this functions passes an invalid ContainerOption.
   * @param vec one-hot-encoded vector
   * @return decoded FeatureVector
   */
  static FeatureVector oneHotDecode(Eigen::VectorXd vec) {
    if (static_cast<size_t>(vec.size()) != oneHotDims) {
      utils::ExceptionHandler::exception("FeatureVector.oneHotDecode: Expected size {}, got {}", oneHotDims,
                                         vec.size());
    }

    size_t pos = 0;
    double cellSizeFactor = vec[pos++];

    // get traversal
    std::optional<TraversalOption> traversal{};
    for (auto &[option, _] : TraversalOption::getOptionNames()) {
      if (vec[pos++] == 1.) {
        if (traversal) {
          utils::ExceptionHandler::exception(
              "FeatureVector.oneHotDecode: Vector encodes more than one traversal. (More than one value for traversal "
              "equals 1.)");
        }
        traversal = option;
      }
    }
    if (not traversal) {
      utils::ExceptionHandler::exception(
          "FeatureVector.oneHotDecode: Vector encodes no traversal. (All values for traversal equal 0.)");
    }

    // get load estimator
    std::optional<LoadEstimatorOption> loadEstimator{};
    for (auto &[option, _] : LoadEstimatorOption::getOptionNames()) {
      if (vec[pos++] == 1.) {
        if (loadEstimator) {
          utils::ExceptionHandler::exception(
              "FeatureVector.oneHotDecode: Vector encodes more than one load estimator. (More than one value for load "
              "estimator "
              "equals 1.)");
        }
        loadEstimator = option;
      }
    }
    if (not loadEstimator) {
      utils::ExceptionHandler::exception(
          "FeatureVector.oneHotDecode: Vector encodes no load estimator. (All values for load estimator equal 0.)");
    }

    // get data layout
    std::optional<DataLayoutOption> dataLayout = {};
    for (auto &[option, _] : DataLayoutOption::getOptionNames()) {
      if (vec[pos++] == 1.) {
        if (dataLayout) {
          utils::ExceptionHandler::exception(
              "FeatureVector.oneHotDecode: Vector encodes more than one data layout. (More than one value for "
              "dataLayout equals 1.)");
        }
        dataLayout = option;
      }
    }
    if (not dataLayout) {
      utils::ExceptionHandler::exception(
          "FeatureVector.oneHotDecode: Vector encodes no data layout. (All values for dataLayout equal 0.)");
    }

    // get newton3
    std::optional<Newton3Option> newton3 = {};
    for (auto &[option, _] : Newton3Option::getOptionNames()) {
      if (vec[pos++] == 1.) {
        if (newton3) {
          utils::ExceptionHandler::exception(
              "FeatureVector.oneHotDecode: Vector encodes more than one newton3. (More than one value for newton3 "
              "equals 1.)");
        }
        newton3 = option;
      }
    }
    if (not newton3) {
      utils::ExceptionHandler::exception(
          "FeatureVector.oneHotDecode: Vector encodes no newton3. (All values for newton3 equal 0.)");
    }

    return FeatureVector(ContainerOption(), cellSizeFactor, *traversal, *loadEstimator, *dataLayout, *newton3);
  }

  /**
   * Create n latin-hypercube-samples from given featureSpace.
   * Container Option of samples are set to -1, because tuning currently
   * ignores this option.
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @param traversals
   * @param loadEstimators
   * @param dataLayouts
   * @param newton3
   * @return vector of sample featureVectors
   */
  static std::vector<FeatureVector> lhsSampleFeatures(size_t n, Random &rng, const NumberSet<double> &cellSizeFactors,
                                                      const std::set<TraversalOption> &traversals,
                                                      const std::set<LoadEstimatorOption> &loadEstimators,
                                                      const std::set<DataLayoutOption> &dataLayouts,
                                                      const std::set<Newton3Option> &newton3) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);
    auto tr = rng.uniformSample(traversals, n);
    auto le = rng.uniformSample(loadEstimators, n);
    auto dl = rng.uniformSample(dataLayouts, n);
    auto n3 = rng.uniformSample(newton3, n);

    std::vector<FeatureVector> result;
    for (size_t i = 0; i < n; ++i) {
      result.emplace_back(ContainerOption(), csf[i], tr[i], le[i], dl[i], n3[i]);
    }

    return result;
  }
};

/**
 * Stream insertion operator.
 * @param os
 * @param featureVector
 * @return
 */
inline std::ostream &operator<<(std::ostream &os, const FeatureVector &featureVector) {
  return os << featureVector.toString();
}

}  // namespace autopas
