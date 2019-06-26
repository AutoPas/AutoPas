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
   * Default constructor. Results in invalid vector.
   */
  FeatureVector() : Configuration() {}

  /**
   * Constructor
   * @param _container
   * @param _traversal
   * @param _dataLayout
   * @param _newton3
   * @param _cellSizeFactor
   */
  FeatureVector(ContainerOption _container, double _cellSizeFactor, TraversalOption _traversal,
                DataLayoutOption _dataLayout, Newton3Option _newton3)
      : Configuration(_container, _cellSizeFactor, _traversal, _dataLayout, _newton3) {}

  /**
   * Construct from Configuration
   * @param conf
   */
  FeatureVector(Configuration conf) : Configuration(conf) {}

  /**
   * Distance between two FeatureVectors
   * @param other
   * @return
   */
  FeatureVector operator-(const FeatureVector &other) const {
    ContainerOption co = ContainerOption((traversal == other.traversal) ? 0 : 1);
    double cfs = (cellSizeFactor - other.cellSizeFactor);
    TraversalOption to = TraversalOption((traversal == other.traversal) ? 0 : 1);
    DataLayoutOption dlo = DataLayoutOption((dataLayout == other.dataLayout) ? 0 : 1);
    Newton3Option n3o = Newton3Option((newton3 == other.newton3) ? 0 : 1);
    return FeatureVector(co, cfs, to, dlo, n3o);
  }

  /**
   * Cast to Eigen::VectorXd ignoring ContainerOption
   * @return
   */
  operator Eigen::VectorXd() const {
    Eigen::VectorXd result(4);
    result << cellSizeFactor, static_cast<double>(traversal), static_cast<double>(dataLayout),
        static_cast<double>(newton3);

    return result;
  }

  /**
   * Create n latin-hypercube-samples from given featureSpace
   * @param n number of samples
   * @param rng
   * @param cellSizeFactors
   * @param traversals
   * @param dataLayouts
   * @param newton3
   * @return
   */
  static std::vector<FeatureVector> lhsSampleFeatures(size_t n, Random &rng, const NumberSet<double> &cellSizeFactors,
                                                      std::set<TraversalOption> traversals,
                                                      std::set<DataLayoutOption> dataLayouts,
                                                      std::set<Newton3Option> newton3) {
    // create n samples from each set
    auto csf = cellSizeFactors.uniformSample(n, rng);
    auto tr = rng.uniformSample(traversals, n);
    auto dl = rng.uniformSample(dataLayouts, n);
    auto n3 = rng.uniformSample(newton3, n);

    std::vector<FeatureVector> result;
    for (unsigned i = 0; i < n; ++i) {
      result.emplace_back(ContainerOption(-1), csf[i], tr[i], dl[i], n3[i]);
    }

    return result;
  }
};
}  // namespace autopas
