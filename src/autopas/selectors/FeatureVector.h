/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

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
   * Distance between two FeatureVectors. Container is ignored.
   * @param other
   * @return array. distance in each dimension
   */
  std::array<double, 4> operator-(const FeatureVector &other) const {
    std::array<double, 4> result;
    result[0] = (cellSizeFactor - other.cellSizeFactor);
    result[1] = ((traversal == other.traversal) ? 0 : 1);
    result[2] = ((dataLayout == other.dataLayout) ? 0 : 1);
    result[3] = ((newton3 == other.newton3) ? 0 : 1);
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
