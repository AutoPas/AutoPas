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
  static constexpr size_t featureSpaceDims = 4;

  /**
   * Dimensions of a one-hot-encoded vector
   * = 1 (cellSizeFactor) + traversals + dataLayouts + newton3
   */
  inline static size_t oneHotDims = 1 + TraversalOption::getOptionNames().size() +
                                    DataLayoutOption::getOptionNames().size() + Newton3Option::getOptionNames().size();

  /**
   * Default constructor. Results in invalid vector.
   */
  FeatureVector() : Configuration() {}

  /**
   * Constructor
   * @param container
   * @param traversal
   * @param dataLayout
   * @param newton3
   * @param cellSizeFactor
   */
  FeatureVector(ContainerOption container, double cellSizeFactor, TraversalOption traversal,
                DataLayoutOption dataLayout, Newton3Option newton3)
      : Configuration(container, cellSizeFactor, traversal, dataLayout, newton3) {}

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
        dataLayout == other.dataLayout ? 0. : 1., newton3 == other.newton3 ? 0. : 1.;
    return result;
  }

  /**
   * Cast to Eigen::VectorXd ignoring ContainerOption.
   * @return
   */
  operator Eigen::VectorXd() const {
    Eigen::VectorXd result(featureSpaceDims);
    result << cellSizeFactor, static_cast<double>(traversal), static_cast<double>(dataLayout),
        static_cast<double>(newton3);

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

    return FeatureVector(ContainerOption(), cellSizeFactor, *traversal, *dataLayout, *newton3);
  }

  /**
   * Convert Feature vector to cluster representation for GaussianCluster.
   * Discrete values are encoded using their index in given std::vector.
   * @param containerTraversalOptions allowed container-traversal pairs
   * @param dataLayoutOptions allowed data layouts
   * @param newton3Options allowed newton3 options
   * @return cluster encoded vector
   */
  [[nodiscard]] std::pair<Eigen::VectorXi, Eigen::VectorXd> convertToCluster(
      const std::vector<std::pair<ContainerOption, TraversalOption>> &containerTraversalOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options) const {
    int containerTraversalIndex = static_cast<int>(std::distance(
        containerTraversalOptions.begin(), std::find(containerTraversalOptions.begin(), containerTraversalOptions.end(),
                                                     std::make_pair(container, traversal))));
    int dataLayoutIndex = static_cast<int>(std::distance(
        dataLayoutOptions.begin(), std::find(dataLayoutOptions.begin(), dataLayoutOptions.end(), dataLayout)));
    int newton3Index = static_cast<int>(
        std::distance(newton3Options.begin(), std::find(newton3Options.begin(), newton3Options.end(), newton3)));

    Eigen::Vector3i vecDiscrete({containerTraversalIndex, dataLayoutIndex, newton3Index});
    Eigen::VectorXd vecContinuous(1);
    vecContinuous << cellSizeFactor;
    return std::make_pair(vecDiscrete, vecContinuous);
  }

  /**
   * Inverse of convertToCluster. Convert cluster representation back
   * to Feature vector.
   * @param vec cluster encoded vector
   * @param containerTraversalOptions allowed container-traversal pairs
   * @param dataLayoutOptions allowed data layouts
   * @param newton3Options allowed newton3 options
   * @return decoded vector
   */
  static FeatureVector convertFromCluster(
      const std::pair<Eigen::VectorXi, Eigen::VectorXd> &vec,
      const std::vector<std::pair<ContainerOption, TraversalOption>> &containerTraversalOptions,
      const std::vector<DataLayoutOption> &dataLayoutOptions, const std::vector<Newton3Option> &newton3Options) {
    const auto &[vecDiscrete, vecContinuous] = vec;
    auto [container, traversal] = containerTraversalOptions[vecDiscrete[0]];

    return FeatureVector(container, vecContinuous[0], traversal, dataLayoutOptions[vecDiscrete[1]],
                         newton3Options[vecDiscrete[2]]);
  }

  /**
   * Get cluster-encoded neighbours of given target.
   * Neighbours are all configurations which differ in at most one configuration from target
   * @param target
   * @param dimRestrictions restriction on each dimension
   * @return all neighbours
   */
  static std::vector<Eigen::VectorXi> neighboursManhattan1(const Eigen::VectorXi &target,
                                                           const std::vector<int> &dimRestrictions) {
    std::vector<Eigen::VectorXi> result;
    // neighbours should contain #(possible values for each dimension) - #dimensions (initial vector is skipped once per
    // dimension)
    result.reserve(std::accumulate(dimRestrictions.begin(), dimRestrictions.end(), -dimRestrictions.size()));

    // for each dimension
    for (int i = 0; i < target.size(); ++i) {
      // initial value
      auto init = target[i];

      // for each possible value of that dimension
      for (int x = 0; x < dimRestrictions[i]; ++x) {
        // skip initial value
        if (x != init) {
          auto neighbour = target;
          neighbour[i] = x;
          result.push_back(std::move(neighbour));
        }
      }
    }
    return result;
  }

  /**
   * Create n latin-hypercube-samples from given featureSpace.
   * Container Option of samples are set to -1, because tuning currently
   * ignores this option.
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @param traversals
   * @param dataLayouts
   * @param newton3
   * @return vector of sample featureVectors
   */
  template <typename TraversalContainer, typename DataLayoutContainer, typename Newton3Container>
  static typename std::enable_if_t<
      std::conjunction_v<std::is_same<typename TraversalContainer::value_type, TraversalOption>,
                         std::is_same<typename DataLayoutContainer::value_type, DataLayoutOption>,
                         std::is_same<typename Newton3Container::value_type, Newton3Option>>,
      std::vector<FeatureVector>>
  lhsSampleFeatures(size_t n, Random &rng, const NumberSet<double> &cellSizeFactors,
                    const TraversalContainer &traversals, const DataLayoutContainer &dataLayouts,
                    const Newton3Container &newton3) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);
    auto tr = rng.uniformSample(traversals.begin(), traversals.end(), n);
    auto dl = rng.uniformSample(dataLayouts.begin(), dataLayouts.end(), n);
    auto n3 = rng.uniformSample(newton3.begin(), newton3.end(), n);

    std::vector<FeatureVector> result;
    for (size_t i = 0; i < n; ++i) {
      result.emplace_back(ContainerOption(), csf[i], tr[i], dl[i], n3[i]);
    }

    return result;
  }
  /**
   * Create n latin-hypercube-samples from given featureSpace only considering continuous values.
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @return vector of sample featureVectors
   */
  static std::vector<Eigen::VectorXd> lhsSampleFeatureContinuous(size_t n, Random &rng,
                                                                 const NumberSet<double> &cellSizeFactors) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);

    std::vector<Eigen::VectorXd> result;
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      Eigen::VectorXd vec(1);
      vec << csf[i];
      result.emplace_back(vec);
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
