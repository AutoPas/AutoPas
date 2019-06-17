/**
 * @file Feature.cpp
 * @author Jan Nguyen
 * @date 22.05.19
 */

#include "Feature.h"

namespace autopas {

Feature::~Feature() = default;

double DoubleFeature::operator-(const Feature &other) const {
  return _value - static_cast<const DoubleFeature &>(other)._value;
}
std::unique_ptr<Feature> DoubleFeature::clone() const { return std::make_unique<DoubleFeature>(_value); }

double EnumFeature::operator-(const Feature &other) const {
  // if same return 0 else 1
  return (_value == static_cast<const EnumFeature &>(other)._value) ? 0. : 1.;
}
std::unique_ptr<Feature> EnumFeature::clone() const { return std::make_unique<EnumFeature>(_value); }

}  // namespace autopas
