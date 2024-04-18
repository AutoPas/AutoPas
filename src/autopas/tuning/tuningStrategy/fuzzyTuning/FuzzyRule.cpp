/**
 * @file FuzzyRule.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "FuzzyRule.h"

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::fuzzy_logic {

FuzzyRule::FuzzyRule(const FuzzySet &antecedent, const FuzzySet &consequent)
    : _antecedent(antecedent), _consequent(consequent) {}

FuzzySet FuzzyRule::apply(const std::map<std::string, double> &data) const {
  const double cut = _antecedent.evaluate_membership(data);

  const std::string newLinguisticTerm = fmt::format("{}↑{:0.2f}", _consequent.getLinguisticTerm(), cut);

  const auto cut_consequent = FuzzySet(
      newLinguisticTerm, [cut, this](auto data) { return std::min(cut, _consequent.evaluate_membership(data)); },
      _consequent.getCrispSet());

  return cut_consequent;
}

}  // namespace autopas::fuzzy_logic
