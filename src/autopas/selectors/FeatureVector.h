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
  static constexpr size_t featureSpaceDims = 5;

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
   * @param verletRebuildFrequency
   */
  FeatureVector(ContainerOption container, double cellSizeFactor, TraversalOption traversal,
                LoadEstimatorOption loadEstimator, DataLayoutOption dataLayout, Newton3Option newton3, int verletRebuildFrequency)
      : Configuration(container, cellSizeFactor, traversal, loadEstimator, dataLayout, newton3, verletRebuildFrequency) {}

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
        dataLayout == other.dataLayout ? 0. : 1., newton3 == other.newton3 ? 0. : 1., verletRebuildFrequency- other.verletRebuildFrequency;
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
