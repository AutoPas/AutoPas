/**
 * @file FuzzySet.cpp
 * @author Manuel Lerchner
 * @date 17.04.24
 */

#include "FuzzySet.h"

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::fuzzy_logic {

FuzzySet::FuzzySet(const std::string &linguisticTerm, const MembershipFunction &membershipFunction)
    : _linguisticTerm(linguisticTerm), _membershipFunction(std::move(membershipFunction)) {}

FuzzySet::FuzzySet(const std::string &linguisticTerm, const MembershipFunction &membershipFunction,
                   const std::shared_ptr<CrispSet> &crispSet)
    : _linguisticTerm(linguisticTerm),
      _membershipFunction(membershipFunction),
      _crispSet(crispSet) {}

FuzzySet::FuzzySet(const std::string &linguisticTerm, const FuzzySet::BaseMembershipFunction &membershipFunction)
    : _linguisticTerm(std::move(linguisticTerm)), _baseMembershipFunction(membershipFunction) {}

double FuzzySet::evaluate_membership(const std::map<std::string, double> &data) const {
  if (_baseMembershipFunction.has_value()) {
    // The current fuzzy set is a base set and has therefore a membership function which can be evaluated with a
    // single value.
    const auto crisp_dimensions = _crispSet->getDimensions();
    if (crisp_dimensions.size() != 1) {
      autopas::utils::ExceptionHandler::exception(
          "FuzzySet::evaluate_membership: FuzzySet is a base set and can therefore only consist of a single "
          "dimension.");
    }
    const auto dimension_name = crisp_dimensions.begin()->first;
    return (*_baseMembershipFunction)(data.at(dimension_name));
  } else {
    // The current fuzzy set is a derived set and has needs to delegate the evaluation recursively to its base sets.
    return _membershipFunction(data);
  }
}

double FuzzySet::centroid(size_t numSamples) const {
  const auto crisp_dimensions = _crispSet->getDimensions();
  if (crisp_dimensions.size() != 1) {
    autopas::utils::ExceptionHandler::exception(
        "FuzzySet::centroid: Can only calculate the centroid of a FuzzySet with a single dimension.");
  }

  const auto [dimensionName, range] = *crisp_dimensions.begin();
  const auto [minBoundary, maxBoundary] = range;

  // Uses the formula centroid_x = sum(x*y) / sum(y) to calculate the centroid of the fuzzy set numerically.
  double numerator = 0;
  double denominator = 0;
  for (double x = minBoundary; x <= maxBoundary; x += (maxBoundary - minBoundary) / numSamples) {
    std::map<std::string, double> data = {{dimensionName, x}};
    const auto membership = evaluate_membership(data);
    numerator += x * membership;
    denominator += membership;
  }

  return denominator == 0 ? 0 : numerator / denominator;
}

FuzzySet FuzzySet::operator&&(const FuzzySet &rhs) const {
  const std::string newLinguisticTerm = fmt::format("({} and {})", _linguisticTerm, rhs._linguisticTerm);
  const auto newCrispSet = *_crispSet * *rhs._crispSet;
  auto newMembershipFunction = [this, rhs](auto data) {
    return std::min(_membershipFunction(data), rhs._membershipFunction(data));
  };
  return FuzzySet(newLinguisticTerm, newMembershipFunction, std::make_shared<CrispSet>(newCrispSet));
}

FuzzySet FuzzySet::operator||(const FuzzySet &rhs) const {
  const std::string newLinguisticTerm = fmt::format("({} or {})", _linguisticTerm, rhs._linguisticTerm);
  const auto newCrispSet = (*_crispSet) * (*rhs._crispSet);
  auto newMembershipFunction = [this, rhs](auto data) {
    return std::max(evaluate_membership(data), rhs.evaluate_membership(data));
  };
  return FuzzySet(newLinguisticTerm, newMembershipFunction, std::make_shared<CrispSet>(newCrispSet));
}

FuzzySet FuzzySet::operator!() const {
  const std::string newLinguisticTerm = fmt::format("not {}", _linguisticTerm);
  const auto newCrispSet = *_crispSet;
  auto newMembershipFunction = [this](auto data) { return 1 - _membershipFunction(data); };
  return FuzzySet(newLinguisticTerm, newMembershipFunction, std::make_shared<CrispSet>(newCrispSet));
}

const std::string &FuzzySet::getLinguisticTerm() const { return _linguisticTerm; }

const FuzzySet::MembershipFunction &FuzzySet::getMembershipFunction() const { return _membershipFunction; }

const std::shared_ptr<CrispSet> &FuzzySet::getCrispSet() const { return _crispSet; }

void FuzzySet::setCrispSet(const std::shared_ptr<CrispSet> &crispSet) { _crispSet = crispSet; }

void FuzzySet::setBaseMembershipFunction(const FuzzySet::BaseMembershipFunction &baseMembershipFunction) {
  _baseMembershipFunction = baseMembershipFunction;
}

}  // namespace autopas::fuzzy_logic
