/**
 * @file FuzzySet.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzySet.h"

#include <iostream>
#include <numeric>
#include <utility>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/Math.h"

namespace autopas::FuzzyLogic {

FuzzySet::FuzzySet(std::string linguisticTerm, BaseMembershipFunction &&baseMembershipFunction)
    : _linguisticTerm(std::move(linguisticTerm)), _membershipFunction(std::move(baseMembershipFunction)) {}

FuzzySet::FuzzySet(std::string linguisticTerm, ComposedMembershipFunction &&membershipFunction,
                   const std::shared_ptr<CrispSet> &crispSet)
    : _linguisticTerm(std::move(linguisticTerm)),
      _membershipFunction(std::move(membershipFunction)),
      _crispSet(crispSet) {}

double FuzzySet::evaluate_membership(const Data &data) const {
  if (std::holds_alternative<BaseMembershipFunction>(_membershipFunction)) {
    // The current fuzzy set is a base set and has therefore a membership function which can be evaluated with a
    // single value.
    const auto crisp_dimensions = _crispSet->getDimensions();
    if (crisp_dimensions.size() != 1) {
      autopas::utils::ExceptionHandler::exception("A base fuzzy set can only have one dimension as input");
    }
    const std::string dimension_name = crisp_dimensions.begin()->first;
    // Extract the value for the current dimension and evaluate the membership function.

    if (data.find(dimension_name) == data.end()) {
      autopas::utils::ExceptionHandler::exception("The data does not contain a value for the dimension " +
                                                  dimension_name);
    }
    return std::get<2>(std::get<BaseMembershipFunction>(_membershipFunction))(data.at(dimension_name));
  } else if (std::holds_alternative<ComposedMembershipFunction>(_membershipFunction)) {
    // The current fuzzy set is a derived set and has needs to delegate the evaluation recursively to its base sets.
    return std::get<ComposedMembershipFunction>(_membershipFunction)(data);
  } else {
    autopas::utils::ExceptionHandler::exception("Unknown membership function type");
    throw std::runtime_error("Unknown membership function type");
  }
}

double FuzzySet::defuzzify(DefuzzificationMethodOption method, size_t numSamples) const {
  switch (method) {
    case DefuzzificationMethodOption::CoG:
      return centerOfGravity(numSamples);
    case DefuzzificationMethodOption::MoM:
      return meanOfMaximum(numSamples);
    default:
      autopas::utils::ExceptionHandler::exception("Unknown defuzzification method");
      throw std::runtime_error("Unknown defuzzification method");
  }
}

double FuzzySet::centerOfGravity(size_t numSamples) const {
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

double FuzzySet::meanOfMaximum(size_t numSamples) const {
  const auto crisp_dimensions = _crispSet->getDimensions();

  if (crisp_dimensions.size() != 1) {
    autopas::utils::ExceptionHandler::exception(
        "The mean of meanOfMaximum can only be calculated for one-dimensional fuzzy sets");
  }

  const auto [dimensionName, range] = *crisp_dimensions.begin();
  const auto [minBoundary, maxBoundary] = range;

  // finds the x-coordinate of the maximum of the fuzzy set numerically.
  double maxMembershipValue = 0;

  std::vector<double> maxMembershipXs;

  std::map<std::string, double> data = {{dimensionName, 0.0}};
  for (double x = minBoundary; x <= maxBoundary; x += (maxBoundary - minBoundary) / (numSamples - 1)) {
    data[dimensionName] = x;
    const double membership = evaluate_membership(data);
    if (membership > maxMembershipValue) {
      maxMembershipValue = membership;
      maxMembershipXs.clear();
    }
    // To avoid missing a maximum due to numerical inaccuracies, we use a small tolerance to compare the two values.
    if (utils::Math::isNear(membership, maxMembershipValue, 1e-5)) {
      maxMembershipXs.push_back(x);
    }
  }

  return maxMembershipXs.empty() ? 0
                                 : std::accumulate(maxMembershipXs.begin(), maxMembershipXs.end(), 0.0) /
                                       static_cast<double>(maxMembershipXs.size());
}

std::string FuzzySet::printBaseMembershipFunction() const {
  auto baseMembershipFunction = std::get<BaseMembershipFunction>(_membershipFunction);
  std::string parameters = std::accumulate(
      std::get<1>(baseMembershipFunction).begin(), std::get<1>(baseMembershipFunction).end(), std::string(""),
      [](const std::string &acc, const double &b) { return acc + std::to_string(b) + ", "; });

  return fmt::format(R"("{}": {}({}))", _linguisticTerm, std::get<0>(baseMembershipFunction),
                     parameters.substr(0, parameters.size() - 2));
}

FuzzySet::operator std::string() const {
  if (std::holds_alternative<BaseMembershipFunction>(_membershipFunction)) {
    // If the fuzzy set is a base set, the membership function is printed in the form "dimension == value".
    const auto [dimensionName, _] = *_crispSet->getDimensions().begin();
    return fmt::format(R"("{}" == "{}")", dimensionName, _linguisticTerm);
  } else {
    // If the fuzzy set is a derived set, just print the linguistic term.
    return _linguisticTerm;
  }
}

const std::string &FuzzySet::getLinguisticTerm() const { return _linguisticTerm; }

const std::shared_ptr<CrispSet> &FuzzySet::getCrispSet() const { return _crispSet; }

void FuzzySet::setCrispSet(const std::shared_ptr<CrispSet> &crispSet) { _crispSet = crispSet; }

std::shared_ptr<FuzzySet> operator||(const std::shared_ptr<FuzzySet> &lhs, const std::shared_ptr<FuzzySet> &rhs) {
  const std::string newLinguisticTerm = fmt::format("({} || {})", std::string(*lhs), std::string(*rhs));
  const auto newCrispSet = (*lhs->_crispSet) * (*rhs->_crispSet);
  auto newMembershipFunction = [lhs, rhs](auto data) {
    return std::max((lhs->evaluate_membership(data)), (rhs->evaluate_membership(data)));
  };
  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(newMembershipFunction), newCrispSet);
}

std::shared_ptr<FuzzySet> operator&&(const std::shared_ptr<FuzzySet> &lhs, const std::shared_ptr<FuzzySet> &rhs) {
  const std::string newLinguisticTerm = fmt::format("({} && {})", std::string(*lhs), std::string(*rhs));
  const auto newCrispSet = (*lhs->_crispSet) * (*rhs->_crispSet);
  auto newMembershipFunction = [lhs, rhs](auto data) {
    return std::min((lhs->evaluate_membership(data)), (rhs->evaluate_membership(data)));
  };
  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(newMembershipFunction), newCrispSet);
}

std::shared_ptr<FuzzySet> operator!(const std::shared_ptr<FuzzySet> &fuzzySet) {
  const std::string newLinguisticTerm = fmt::format("!({})", std::string(*fuzzySet));
  const auto newCrispSet = fuzzySet->_crispSet;
  auto newMembershipFunction = [fuzzySet](auto data) { return 1 - (fuzzySet->evaluate_membership(data)); };
  return std::make_shared<FuzzySet>(newLinguisticTerm, std::move(newMembershipFunction), newCrispSet);
}
}  // namespace autopas::FuzzyLogic
