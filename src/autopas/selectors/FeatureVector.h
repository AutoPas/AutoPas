/**
 * @file FeatureVector.h
 * @author Jan Nguyen
 * @date 17.05.19
 */

#pragma once

#include <Eigen/Dense>
#include <memory>
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * Vector containing any number of features
 */
class FeatureVector {
  /**
   * Class describing a feature
   */
  class Feature {
   public:
    virtual ~Feature() = default;

    /**
     * Subtract two features
     * @param other
     * @return feature distance
     */
    virtual double operator-(const Feature&) const;
  };

  /**
   * Feature described by a double
   */
  class DoubleFeature : public Feature {
   private:
    double _value;

   public:
    /**
     * Construct double feature from double
     * @param value
     */
    DoubleFeature(double value) : _value(value) {}

    double operator-(const Feature& other) const override;
  };

 private:
  std::vector<std::shared_ptr<Feature>> _vector;

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

  void addFeature(double doubleFeature) { _vector.push_back(std::make_shared<DoubleFeature>(doubleFeature)); }
};

double FeatureVector::Feature::operator-(const Feature&) const { return 0.; }
double FeatureVector::DoubleFeature::operator-(const Feature& other) const {
  return _value - dynamic_cast<const DoubleFeature&>(other)._value;
}
}  // namespace autopas
