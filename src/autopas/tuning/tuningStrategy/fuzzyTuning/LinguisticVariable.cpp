/**
 * @file LinguisticVariable.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "LinguisticVariable.h"

namespace autopas::fuzzy_logic {

LinguisticVariable::LinguisticVariable(const std::string &name, const std::pair<double, double> &range) {
  _crispSet = std::make_shared<CrispSet>(name, range);
}

void LinguisticVariable::addLinguisticTerm(const std::shared_ptr<FuzzySet> &fuzzySet) {
  fuzzySet->setCrispSet(_crispSet);
  _linguisticTerms[fuzzySet->getLinguisticTerm()] = fuzzySet;
}

std::shared_ptr<FuzzySet> LinguisticVariable::operator==(const std::string &linguisticTerm) const {
  return _linguisticTerms.at(linguisticTerm);
}

}  // namespace autopas::fuzzy_logic
