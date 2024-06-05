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

struct FuzzyControlSettings {
  /**
   * The defuzzification method to use. Default is "centroid".
   */
  std::string defuzzificationMethod = "centroid";
};

class FuzzyControlSystem {
 public:
  /**
   * Constructs an empty FuzzyControlSystem with the given settings.
   */
  explicit FuzzyControlSystem(FuzzyControlSettings settings = {});

  /**
   * Adds a new FuzzyRule to the FuzzyControlSystem.
   * @param rule The FuzzyRule to add.
   *
   * Additionally, the method checks if all the consequents of the FuzzyControlSystem have the same domain.
   * Otherwise, an exception is thrown.
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

  /**
   * Returns a string representation of the FuzzyControlSystem.
   */
  explicit operator std::string() const;

 private:
  /**
   * The settings of the FuzzyControlSystem.
   */
  FuzzyControlSettings _settings;

  /**
   * All rules of the FuzzyControlSystem.
   */
  std::vector<FuzzyRule> _rules;

  /**
   * The output domain of the FuzzyControlSystem.
   */
  std::optional<std::string> _outputDomain;
};

}  // namespace autopas::fuzzy_logic