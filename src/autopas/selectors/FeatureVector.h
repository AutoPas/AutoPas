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
 * FeatureVector representation of a Configuration.
 */
class FeatureVector : public Configuration {
 public:
  /**
   * Number of tune-able dimensions.
   * container-traversal-estimator + dataLayout + newton3 + cellSizeFactor
   */
  static constexpr size_t featureSpaceDims = 4;

  /**
   * Number of tune-able continuous dimensions.
   */
  static constexpr size_t featureSpaceContinuousDims = 1;

  /**
   * Consider Container, Traversal and LoadEstimator options as one dimension.
   */
  using ContainerTraversalEstimatorOption = std::tuple<ContainerOption, TraversalOption, LoadEstimatorOption>;

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
   * @param other
   * @return
   */
  Eigen::VectorXd operator-(const FeatureVector &other) const {
    Eigen::VectorXd result(featureSpaceDims);
    result << cellSizeFactor - other.cellSizeFactor,
        (container == other.container and traversal == other.traversal and loadEstimator == other.loadEstimator) ? 0.
                                                                                                                 : 1.,
        dataLayout == other.dataLayout ? 0. : 1., newton3 == other.newton3 ? 0. : 1.;
    return result;
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
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @param containerTraversalEstimators
   * @param dataLayouts
   * @param newton3
   * @return vector of sample featureVectors
   */
  template <typename ContainerTraversalEstimatorContainer, typename DataLayoutContainer, typename Newton3Container>
  static typename std::enable_if_t<
      std::conjunction_v<
          std::is_same<typename ContainerTraversalEstimatorContainer::value_type, ContainerTraversalEstimatorOption>,
          std::is_same<typename DataLayoutContainer::value_type, DataLayoutOption>,
          std::is_same<typename Newton3Container::value_type, Newton3Option>>,
      std::vector<FeatureVector>>
  lhsSampleFeatures(size_t n, Random &rng, const NumberSet<double> &cellSizeFactors,
                    const ContainerTraversalEstimatorContainer &containerTraversalEstimators,
                    const DataLayoutContainer &dataLayouts, const Newton3Container &newton3) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);
    auto cte = rng.uniformSample(containerTraversalEstimators.begin(), containerTraversalEstimators.end(), n);
    auto dl = rng.uniformSample(dataLayouts.begin(), dataLayouts.end(), n);
    auto n3 = rng.uniformSample(newton3.begin(), newton3.end(), n);

    std::vector<FeatureVector> result;
    for (size_t i = 0; i < n; ++i) {
      const auto &[container, traversal, estimator] = cte[i];
      result.emplace_back(container, csf[i], traversal, estimator, dl[i], n3[i]);
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
      Eigen::VectorXd vec(featureSpaceContinuousDims);
      vec << csf[i];
      result.emplace_back(vec);
    }

    return result;
  }

  /**
   * Create n latin-hypercube-samples from given featureSpace only considering continuous values and append a value
   * representing the current iteration.
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @param iteration current iteration which may be scaled by some factor
   * @return vector of sample featureVectors
   */
  static std::vector<Eigen::VectorXd> lhsSampleFeatureContinuousWithIteration(size_t n, Random &rng,
                                                                              const NumberSet<double> &cellSizeFactors,
                                                                              double iteration) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);

    std::vector<Eigen::VectorXd> result;
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      Eigen::VectorXd vec(featureSpaceContinuousDims + 1);
      vec << csf[i], iteration;
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
