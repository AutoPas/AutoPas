/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <memory>
#include <vector>
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/Random.h"

namespace autopas {

/**
 * Vector containing any number of features
 */
class FeatureVector : public Configuration {
 public:
  /**
   * Default constructor
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
   * Constructor from configuration
   * @param conf
   */
  FeatureVector(Configuration conf) : Configuration(conf) {}

  /**
   * Distance vector between two feature vectors
   * @param other
   * @return vector<double> distance in each dimension
   */
  std::vector<double> operator-(const FeatureVector &other) const {
    std::vector<double> result;
    result.push_back(cellSizeFactor - other.cellSizeFactor);
    result.push_back((traversal == other.traversal) ? 0 : 1);
    result.push_back((dataLayout == other.dataLayout) ? 0 : 1);
    result.push_back((newton3 == other.newton3) ? 0 : 1);
    return result;
  }

  /**
   * Create equidistant values in given range and
   * append them randomly to given vectors.
   * @param vectors
   * @param featureSet
   * @param rng random number generator
   */
  static void lhsSetCellSizeFactors(std::vector<FeatureVector> &vectors, const NumberSet<double> &featureSet,
                                    Random &rng) {
    // create n samples from given set
    auto pool = featureSet.uniformSample(vectors.size(), rng);

    // set the features
    for (unsigned i = 0; i < vectors.size(); ++i) {
      vectors[i].cellSizeFactor = pool[i];
    }
  }

  /**
   * Randomly set the TraversalOption of all FeatureVectors with
   * random values from the featureSpace.
   *
   * @param vectors
   * @param featureSpace
   * @param rng random number generator
   */
  static void lhsSetTraversals(std::vector<FeatureVector> &vectors, std::set<TraversalOption> featureSpace,
                               Random &rng) {
    // create n samples from the feature space
    auto pool = rng.uniformSample(featureSpace, vectors.size());

    // set the features
    for (unsigned i = 0; i < vectors.size(); ++i) {
      vectors[i].traversal = pool[i];
    }
  }

  /**
   * Randomly set the DataLayoutOption of all FeatureVectors with
   * random values from the featureSpace.
   *
   * @param vectors
   * @param featureSpace
   * @param rng random number generator
   */
  static void lhsSetDataLayouts(std::vector<FeatureVector> &vectors, std::set<DataLayoutOption> featureSpace,
                                Random &rng) {
    // create n samples from the feature space
    auto pool = rng.uniformSample(featureSpace, vectors.size());

    // set the features
    for (unsigned i = 0; i < vectors.size(); ++i) {
      vectors[i].dataLayout = pool[i];
    }
  }

  /**
   * Randomly set the Newton3Option of all FeatureVectors with
   * random values from the featureSpace.
   *
   * @param vectors
   * @param featureSpace
   * @param rng random number generator
   */
  static void lhsSetNewton3(std::vector<FeatureVector> &vectors, std::set<Newton3Option> featureSpace, Random &rng) {
    // create n samples from the feature space
    auto pool = rng.uniformSample(featureSpace, vectors.size());

    // set the features
    for (unsigned i = 0; i < vectors.size(); ++i) {
      vectors[i].newton3 = pool[i];
    }
  }
};
}  // namespace autopas
