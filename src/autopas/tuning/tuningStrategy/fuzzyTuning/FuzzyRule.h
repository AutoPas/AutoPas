/**
 * @file FuzzyRule.h
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#pragma once

#include <memory>
#include <vector>

#include "FuzzySet.h"

namespace autopas::fuzzy_logic {

class FuzzyRule {
 public:
  /**
   * Constructs a FuzzyRule of the form: IF antecedent THEN consequent.
   * @param antecedent A FuzzySet representing the antecedent of the FuzzyRule.
   * @param consequent A FuzzySet representing the consequent of the FuzzyRule.
   */
  FuzzyRule(const std::shared_ptr<FuzzySet> &antecedent, const std::shared_ptr<FuzzySet> &consequent);

  /**
   * Applies the FuzzyRule to the given data.
   * @param data A map of the form {dimension_name: value}.
   * @return The cut FuzzySet resulting from the application of the FuzzyRule to the given data.
   */
  [[nodiscard]] std::shared_ptr<FuzzySet> apply(const std::map<std::string, double> &data) const;

 private:
  /**
   * The antecedent of the FuzzyRule.
   */
  const std::shared_ptr<FuzzySet> _antecedent;

  /**
   * The consequent of the FuzzyRule.
   */
  const std::shared_ptr<FuzzySet> _consequent;
};

}  // namespace autopas::fuzzy_logic