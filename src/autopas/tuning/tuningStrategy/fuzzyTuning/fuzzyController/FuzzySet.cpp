/**
 * @file FuzzySet.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzySet.h"

#include <iostream>
#include <utility>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::fuzzy_logic {

FuzzySet::FuzzySet(std::string linguisticTerm, BaseMembershipFunction &&baseMembershipFunction)
    : _linguisticTerm(std::move(linguisticTerm)), _baseMembershipFunction(std::move(baseMembershipFunction)) {}

FuzzySet::FuzzySet(std::string linguisticTerm, ComposedMembershipFunction &&membershipFunction,
                   const std::shared_ptr<CrispSet> &crispSet)
    : _linguisticTerm(std::move(linguisticTerm)),
      _membershipFunction(std::move(membershipFunction)),
      _crispSet(crispSet) {}

double FuzzySet::evaluate_membership(const Data &data) const {
  if (_baseMembershipFunction.has_value()) {
    // The current fuzzy set is a base set and has therefore a membership function which can be evaluated with a
    // single value.
    const auto crisp_dimensions = _crispSet->getDimensions();
    if (crisp_dimensions.size() != 1) {
      autopas::utils::ExceptionHandler::exception("A base fuzzy set can only have one dimension as input");
    }
    const std::string dimension_name = crisp_dimensions.begin()->first;
    // Extract the value for the current dimension and evaluate the membership function.
    return (*_baseMembershipFunction)(data.at(dimension_name));
  } else {
    // The current fuzzy set is a derived set and has needs to delegate the evaluation recursively to its base sets.
    return _membershipFunction(data);
  }
}

double FuzzySet::centroid(size_t numSamples) const {
  const auto crisp_dimensions = _crispSet->getDimensions();

  if (crisp_dimensions.size() != 1) {
    autopas::utils::ExceptionHandler::exception("The centroid can only be calculated for one-dimensional fuzzy sets");
  }

  const auto [dimensionName, range] = *crisp_dimensions.begin();
  const auto [minBoundary, maxBoundary] = range;

  // Uses the formula centroid_x = sum(x*y) / sum(y) to calculate the centroid of the fuzzy set numerically.
  double numerator = 0;
  double denominator = 0;
  std::map<std::string, double> data = {{dimensionName, 0.0}};
  for (double x = minBoundary; x <= maxBoundary; x += (maxBoundary - minBoundary) / (numSamples - 1)) {
    data[dimensionName] = x;
    const double membership = evaluate_membership(data);
    numerator += x * membership;
    denominator += membership;
  }

  return denominator == 0 ? 0 : numerator / denominator;
}

const std::string &FuzzySet::getLinguisticTerm() const { return _linguisticTerm; }

const std::shared_ptr<CrispSet> &FuzzySet::getCrispSet() const { return _crispSet; }

void FuzzySet::setCrispSet(const std::shared_ptr<CrispSet> &crispSet) { _crispSet = crispSet; }

std::shared_ptr<FuzzySet> operator||(const std::shared_ptr<FuzzySet> &lhs, const std::shared_ptr<FuzzySet> &rhs) {
  const std::string newLinguisticTerm = fmt::format("({} || {})", lhs->_linguisticTerm, rhs->_linguisticTerm);
  const auto newCrispSet = (*lhs->_crispSet) * (*rhs->_crispSet);
  auto newMembershipFunction = [lhs, rhs](auto data) {
    return std::max((lhs->evaluate_membership(data)), (rhs->evaluate_membership(data)));
  };
  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(newMembershipFunction), newCrispSet);
}

std::shared_ptr<FuzzySet> operator&&(const std::shared_ptr<FuzzySet> &lhs, const std::shared_ptr<FuzzySet> &rhs) {
  const std::string newLinguisticTerm = fmt::format("({} && {})", lhs->_linguisticTerm, rhs->_linguisticTerm);
  const auto newCrispSet = (*lhs->_crispSet) * (*rhs->_crispSet);
  auto newMembershipFunction = [lhs, rhs](auto data) {
    return std::min((lhs->evaluate_membership(data)), (rhs->evaluate_membership(data)));
  };
  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(newMembershipFunction), newCrispSet);
}

std::shared_ptr<FuzzySet> operator!(const std::shared_ptr<FuzzySet> &set) {
  const std::string newLinguisticTerm = fmt::format("!{}", set->_linguisticTerm);
  const auto newCrispSet = set->_crispSet;
  auto newMembershipFunction = [set](auto data) { return 1 - (set->evaluate_membership(data)); };
  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(newMembershipFunction), newCrispSet);
}
}  // namespace autopas::fuzzy_logic
