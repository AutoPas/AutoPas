/**
 * @file FuzzyRule.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "FuzzyRule.h"

namespace autopas::fuzzy_logic {

FuzzyRule::FuzzyRule(const std::shared_ptr<FuzzySet> &antecedent, const std::shared_ptr<FuzzySet> &consequent)
    : _antecedent(antecedent), _consequent(consequent) {}

std::shared_ptr<FuzzySet> FuzzyRule::apply(const FuzzySet::Data &data) const {
  // calculate how much the antecedent is fulfilled
  const double cut = _antecedent->evaluate_membership(data);

  const std::string newLinguisticTerm = _consequent->getLinguisticTerm() + "↑" + std::to_string(cut);

  // The cut consequent is the minimum of the calculated cut and the consequent membership function.
  const auto cutConsequent = [cut, this](const FuzzySet::Data &data) {
    return std::min(cut, _consequent->evaluate_membership(data));
  };

  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(cutConsequent), _consequent->getCrispSet());
}

}  // namespace autopas::fuzzy_logic
