/**
 * @file Feature.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <memory>

namespace autopas {

/**
 * Class describing a feature
 */
class Feature {
 public:
  virtual ~Feature();

  /**
   * Subtract two features
   * @param other
   * @return feature distance
   */
  virtual double operator-(const Feature& other) const = 0;

  /**
   * Create a copy of a feature
   * @return
   */
  virtual std::unique_ptr<Feature> clone() const = 0;
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
  std::unique_ptr<Feature> clone() const override;
};

Feature::~Feature() = default;

double DoubleFeature::operator-(const Feature& other) const {
  return _value - dynamic_cast<const DoubleFeature&>(other)._value;
}
std::unique_ptr<Feature> DoubleFeature::clone() const { return std::make_unique<DoubleFeature>(_value); }
}  // namespace autopas
