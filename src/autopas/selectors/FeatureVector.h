/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <Eigen/Dense>
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
  static constexpr size_t _featureSpaceDims = 4;

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
   * @param other
   * @return
   */
  Eigen::VectorXd operator-(const FeatureVector &other) const {
    Eigen::VectorXd result(_featureSpaceDims);
    result << cellSizeFactor - other.cellSizeFactor, traversal == other.traversal ? 0. : 1.,
        dataLayout == other.dataLayout ? 0. : 1., newton3 == other.newton3 ? 0. : 1.;

    return result;
  }

  /**
   * Cast to Eigen::VectorXd ignoring ContainerOption.
   * @return
   */
  operator Eigen::VectorXd() const {
    Eigen::VectorXd result(_featureSpaceDims);
    result << cellSizeFactor, static_cast<double>(traversal), static_cast<double>(dataLayout),
        static_cast<double>(newton3);

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
  static std::vector<FeatureVector> lhsSampleFeatures(size_t n, Random &rng, const NumberSet<double> &cellSizeFactors,
                                                      const std::set<TraversalOption> &traversals,
                                                      const std::set<DataLayoutOption> &dataLayouts,
                                                      const std::set<Newton3Option> &newton3) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);
    auto tr = rng.uniformSample(traversals, n);
    auto dl = rng.uniformSample(dataLayouts, n);
    auto n3 = rng.uniformSample(newton3, n);

    std::vector<FeatureVector> result;
    for (size_t i = 0; i < n; ++i) {
      result.emplace_back(ContainerOption(), csf[i], tr[i], dl[i], n3[i]);
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
