/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <memory>
#include <vector>
#include "autopas/selectors/Feature.h"

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
   * Add a DoubleFeature with given value to the vector
   * @param doubleFeature
   */
  void addFeature(double doubleFeature) { _vector.push_back(std::make_unique<DoubleFeature>(doubleFeature)); }
};
}  // namespace autopas
