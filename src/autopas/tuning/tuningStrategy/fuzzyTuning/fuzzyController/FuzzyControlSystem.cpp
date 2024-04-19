/**
 * @file FuzzyControlSystem.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "FuzzyControlSystem.h"

#include <numeric>

namespace autopas::fuzzy_logic {

void FuzzyControlSystem::addRule(const FuzzyRule &rule) { _rules.push_back(rule); };

std::shared_ptr<FuzzySet> FuzzyControlSystem::applyRules(const FuzzySet::Data &data) const {
  // Applies all the rules of the system to get individual fuzzy sets corresponding to the cut-consequents.
  // Then combines them using the Fuzzy-OR operator to get the union set.
  return std::transform_reduce(_rules.begin() + 1, _rules.end(), _rules[0].apply(data), operator||,
                               [&data](const FuzzyRule &rule) { return rule.apply(data); });
}

double FuzzyControlSystem::predict(const FuzzySet::Data &data, size_t numSamples) const {
  auto unionSet = applyRules(data);
  return unionSet->centroid(numSamples);
}

}  // namespace autopas::fuzzy_logic