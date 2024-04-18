/**
 * @file LinguisticVariable.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "LinguisticVariable.h"

namespace autopas::fuzzy_logic {

LinguisticVariable::LinguisticVariable(const std::shared_ptr<CrispSet> &crispSet) : _crispSet(crispSet) {}

void LinguisticVariable::addLinguisticTerm(const std::shared_ptr<FuzzySet> &fuzzySet) {
  fuzzySet->setCrispSet(_crispSet);
  _linguisticTerms[fuzzySet->getLinguisticTerm()] = fuzzySet;
}

const FuzzySet &LinguisticVariable::operator==(const std::string &linguisticTerm) const {
  return *_linguisticTerms.at(linguisticTerm);
}

}  // namespace autopas::fuzzy_logic
