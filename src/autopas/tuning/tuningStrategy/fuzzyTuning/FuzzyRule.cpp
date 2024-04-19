/**
 * @file FuzzyRule.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "FuzzyRule.h"

namespace autopas::fuzzy_logic {

FuzzyRule::FuzzyRule(const std::shared_ptr<FuzzySet> &antecedent, const std::shared_ptr<FuzzySet> &consequent)
    : _antecedent(antecedent), _consequent(consequent) {}

std::shared_ptr<FuzzySet> FuzzyRule::apply(const std::map<std::string, double> &data) const {
  const double cut = _antecedent->evaluate_membership(data);

  const std::string newLinguisticTerm = _consequent->getLinguisticTerm() + "↑" + std::to_string(cut);

  const auto cons =
      std::make_shared<FuzzySet::ComposedMembershipFunction>([cut, this](const std::map<std::string, double> &data) {
        return std::min(cut, _consequent->evaluate_membership(data));
      });

  return std::make_shared<FuzzySet>(newLinguisticTerm, cons, _consequent->getCrispSet());
}

}  // namespace autopas::fuzzy_logic
