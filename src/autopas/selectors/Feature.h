/**
 * @file Feature.h
 * @author Jan Nguyen
 * @date 22.05.19
 */

#pragma once

#include <memory>
#include <set>
#include <vector>

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

  /**
   * Get double representing the Feature
   * @return
   */
  double getValue() { return _value; }
};

/**
 * Feature described by a enum.
 * Different enums always have distance 1 from each other.
 */
class EnumFeature : public Feature {
 private:
  int _value;

 public:
  /**
   * Construct enum feature from int represntation of enum.
   * @param value
   */
  EnumFeature(int value) : _value(value) {}

  double operator-(const Feature& other) const override;
  std::unique_ptr<Feature> clone() const override;

  /**
   * Get int representing the Feature
   * @return
   */
  int getValue() { return _value; }

  /**
   * Convert a set of enums to a vector of EnumFeature
   */
  template <class enumType>
  static std::vector<EnumFeature> set2Vector(std::set<enumType> enums) {
    std::vector<EnumFeature> result;
    for (auto& e : enums) {
      result.push_back(EnumFeature(e));
    }

    return result;
  }
};
}  // namespace autopas
