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
   * Applies the FuzzyRule using the given data.
   * @param data A map of the form {dimension_name: value}.
   * @return The cut consequent FuzzySet resulting from the application of the FuzzyRule.
   */
  [[nodiscard]] std::shared_ptr<FuzzySet> apply(const FuzzySet::Data &data) const;

  /**
   * Returns the antecedent of the FuzzyRule.
   * @return The antecedent of the FuzzyRule.
   */
  [[nodiscard]] const std::shared_ptr<FuzzySet> &getAntecedent() const;

  /**
   * Returns the consequent of the FuzzyRule.
   * @return The consequent of the FuzzyRule.
   */
  [[nodiscard]] const std::shared_ptr<FuzzySet> &getConsequent() const;

  /**
   * Returns a string representation of the FuzzyRule.
   */
  explicit operator std::string() const;

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