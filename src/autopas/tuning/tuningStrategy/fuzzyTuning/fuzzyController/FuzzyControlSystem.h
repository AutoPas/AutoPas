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

namespace autopas::FuzzyLogic {

/**
 * The settings of a FuzzyControlSystem are a map of key-value pairs. The key is the name of the setting and the value
 * is the value of the setting.
 */
using FuzzyControlSettings = std::map<std::string, std::string>;

/**
 * Used to represent a Fuzzy Control System. A Fuzzy Control System is a collection of FuzzyRules that can be applied to
 * a given input to predict an output.
 */
class FuzzyControlSystem {
 public:
  /**
   * Constructs an empty FuzzyControlSystem with the given settings.
   * @param settings The settings of the FuzzyControlSystem.
   */
  explicit FuzzyControlSystem(std::shared_ptr<FuzzyControlSettings> settings);

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
   * @param numSamples The number of samples to use for the numerical defuzzification. Default is 1000.
   * @return The predicted output of the FuzzyControlSystem for the given data.
   *
   * This method performs the full prediction process of the FuzzyControlSystem. It first applies all the rules of the
   * FuzzyControlSystem to the data and unites all the cut-consequents. The resulting FuzzySet is then defuzzified and
   * the calculated value is returned.
   *
   * The method reads the "defuzzificationMethod" and "numSamples" settings from the provided settings and performs the
   * defuzzification accordingly.
   */
  [[nodiscard]] double predict(const FuzzySet::Data &data, size_t numSamples = 1000) const;

  /**
   * Returns a string representation of the FuzzyControlSystem.
   */
  explicit operator std::string() const;

 private:
  /**
   * The settings of the FuzzyControlSystem.
   */
  std::shared_ptr<FuzzyControlSettings> _settings;

  /**
   * All rules of the FuzzyControlSystem.
   */
  std::vector<FuzzyRule> _rules;

  /**
   * The output domain of the FuzzyControlSystem.
   */
  std::optional<std::string> _outputDomain;
};

}  // namespace autopas::FuzzyLogic