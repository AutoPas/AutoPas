/**
 * @file LinguisticVariable.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "LinguisticVariable.h"

#include <numeric>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::FuzzyLogic {

LinguisticVariable::LinguisticVariable(const std::string &name, std::pair<double, double> range) {
  _name = name;
  _crispSet = std::make_shared<CrispSet>(name, range);
}

void LinguisticVariable::addLinguisticTerm(const std::shared_ptr<FuzzySet> &linguisticTerm) {
  linguisticTerm->setCrispSet(_crispSet);
  _linguisticTerms[linguisticTerm->getLinguisticTerm()] = linguisticTerm;
}

std::shared_ptr<FuzzySet> LinguisticVariable::operator==(const std::string &linguisticTerm) const {
  if (_linguisticTerms.find(linguisticTerm) == _linguisticTerms.end()) {
    autopas::utils::ExceptionHandler::exception("Linguistic term " + linguisticTerm + " not found in variable " +
                                                _name);
  }

  return _linguisticTerms.at(linguisticTerm);
}

const std::string &LinguisticVariable::getName() const { return _name; }

LinguisticVariable::operator std::string() const {
  auto dimensions = _crispSet->getDimensions();
  if (dimensions.size() != 1) {
    autopas::utils::ExceptionHandler::exception(
        "A linguistic variable can only be based on a one-dimensional CrispSet");
  }

  const auto [dimensionName, range] = *dimensions.begin();

  std::string linguisticTermsStr =
      std::accumulate(_linguisticTerms.begin(), _linguisticTerms.end(), std::string(""),
                      [](const std::string &acc, const std::pair<const std::string, std::shared_ptr<FuzzySet>> &b) {
                        return acc + "\t" + b.second->printBaseMembershipFunction() + "\n";
                      });

  return fmt::format("LinguisticVariable: domain: \"{}\" range: ({}, {})\n{}", _name, range.first, range.second,
                     linguisticTermsStr);
}

}  // namespace autopas::FuzzyLogic
