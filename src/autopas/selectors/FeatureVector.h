/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <memory>
#include <random>
#include <vector>
#include "autopas/selectors/Feature.h"
#include "autopas/utils/DoubleSet.h"

namespace autopas {

/**
 * Vector containing any number of features
 */
class FeatureVector {
 private:
  std::vector<std::unique_ptr<Feature>> _vector;

 public:
  /**
   * Default constructor
   */
  FeatureVector() : _vector() {}

  /**
   * Initialize vector from list of doubles
   * @param init
   */
  FeatureVector(std::initializer_list<double> init) : _vector() {
    for (double i : init) {
      addFeature(i);
    }
  }

  /**
   * Copy constructor
   * @param other
   */
  FeatureVector(const FeatureVector& other) : _vector() {
    _vector.reserve(other._vector.size());
    for (auto& f : other._vector) {
      _vector.push_back(f->clone());
    }
  }

  /**
   * Copy assignment
   * @param other
   * @return
   */
  FeatureVector& operator=(const FeatureVector& other) {
    _vector.clear();
    _vector.reserve(other._vector.size());
    for (auto& f : other._vector) {
      _vector.push_back(f->clone());
    }

    return *this;
  }

  /**
   * Distance vector between two feature vectors
   * @param other
   * @return vector<double> distance in each dimension
   */
  std::vector<double> operator-(const FeatureVector& other) const {
    std::vector<double> result;
    for (unsigned i = 0; i < _vector.size(); ++i) {
      result.push_back((*_vector[i]) - (*other._vector[i]));
    }

    return result;
  }

  /**
   * get the feature of given index
   * @param index
   * @return
   */
  Feature& operator[](size_t index) { return *_vector[index]; }

  /**
   * Add a feature to the vector
   * @param feature
   */
  void addFeature(const Feature& feature) { _vector.push_back(feature.clone()); }

  /**
   * Add a DoubleFeature with given value to the vector
   * @param doubleFeature
   */
  void addFeature(double doubleFeature) { _vector.push_back(std::make_unique<DoubleFeature>(doubleFeature)); }

  /**
   * Create equidistant values in given range and
   * append them randomly to given vectors.
   * @param vectors
   * @param featureRange
   * @param seed seed for rng
   */
  static void lhsAddFeature(std::vector<FeatureVector>& vectors, const DoubleSet& featureSet, std::default_random_engine& rng) {
    // create n samples from given set
    auto pool = featureSet.uniformSample(vectors.size(), rng);

    // append to feature vectors
    for (unsigned i = 0; i < vectors.size(); ++i) {
      vectors[i].addFeature(pool[i]);
    }
  }

  /**
   * Create a sample from given feature space and
   * append them to the vectors randomly.
   *
   * @param vectors
   * @param featureSpace
   */
  static void lhsAddFeature(std::vector<FeatureVector>& vectors, std::vector<std::unique_ptr<Feature>> featureSpace,
                            std::default_random_engine& rng) {
    // create n values from given pool
    std::vector<std::unique_ptr<Feature>> pool;
    pool.reserve(vectors.size());

    // first fill the pool with copies of the whole feature space
    unsigned minCopies = vectors.size() / featureSpace.size();
    for (unsigned i = 0; i < minCopies; ++i) {
      for (auto& feature : featureSpace) {
        pool.push_back(feature->clone());
      }
    }
    // fill the rest with random samples
    while (pool.size() < vectors.size()) {
      // get a random index
      std::uniform_int_distribution<std::mt19937::result_type> uniDist(0, vectors.size() - pool.size() - 1);
      auto index = uniDist(rng);

      // move from featureSpace vector to pool vector
      pool.push_back(std::unique_ptr<Feature>(featureSpace[index].release()));
      featureSpace.erase(featureSpace.begin() + index);
    }

    // randomize pool
    std::shuffle(std::begin(pool), std::end(pool), rng);

    // append to feature vectors
    for (unsigned i = 0; i < vectors.size(); ++i) {
      vectors[i].addFeature(*pool[i]);
    }
  }
};
}  // namespace autopas
