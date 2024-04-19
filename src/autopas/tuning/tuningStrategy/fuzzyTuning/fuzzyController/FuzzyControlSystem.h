/**
 * @file FuzzyControlSystem.h
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#pragma once

#include <memory>
#include <vector>

#include "FuzzyRule.h"
#include "FuzzySet.h"

namespace autopas::fuzzy_logic {

class FuzzyControlSystem {
 public:
  /**
   * Constructs an empty FuzzyControlSystem.
   */
  FuzzyControlSystem() = default;

  /**
   * Adds a new FuzzyRule to the FuzzyControlSystem.
   * @param rule The FuzzyRule to add.
   */
  void addRule(const FuzzyRule &rule);

  /**
   * Applies all the rules of the FuzzyControlSystem with the given data and calculates the union of all the
   * cut-consequents.
   * @param data A map of the form {dimension_name: value}.
   * @return The FuzzySet resulting from the application of the FuzzyControlSystem to the given data.
   */
  [[nodiscard]] std::shared_ptr<FuzzySet> applyRules(const FuzzySet::Data &data) const;

  /**
   * Predicts the output of the FuzzyControlSystem for the given data.
   * @param data A map of the form {dimension_name: value}.
   * @param numSamples The number of samples to use for the numerical centroid calculation. Default is 100.
   * @return The predicted output of the FuzzyControlSystem for the given data.
   */
  [[nodiscard]] double predict(const FuzzySet::Data &data, size_t numSamples = 100) const;

 private:
  /**
   * All rules of the FuzzyControlSystem.
   */
  std::vector<FuzzyRule> _rules;
};

}  // namespace autopas::fuzzy_logic