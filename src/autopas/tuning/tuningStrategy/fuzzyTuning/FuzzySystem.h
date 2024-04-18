/**
 * @file FuzzySystem.h
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#pragma once

#include <memory>
#include <vector>

#include "FuzzyRule.h"
#include "FuzzySet.h"

namespace autopas::fuzzy_logic {

class FuzzySystem {
 public:
  /**
   * Adds a new FuzzyRule to the FuzzySystem.
   * @param rule The FuzzyRule to add.
   */
  void addRule(const FuzzyRule &rule);

  /**
   * Applies the FuzzySystem to the given data.
   * @param data A map of the form {dimension_name: value}.
   * @return The FuzzySet resulting from the application of the FuzzySystem to the given data.
   */
  FuzzySet applyRules(const std::map<std::string, double> &data) const;

  /**
   * Predicts the output of the FuzzySystem for the given data.
   * @param data A map of the form {dimension_name: value}.
   * @param numSamples The number of samples to use for the centroid calculation.
   * @return The predicted output of the FuzzySystem for the given data.
   */
  double predict(const std::map<std::string, double> &data, size_t numSamples = 100) const;

 private:
  std::vector<FuzzyRule> _rules;
};

}  // namespace autopas::fuzzy_logic