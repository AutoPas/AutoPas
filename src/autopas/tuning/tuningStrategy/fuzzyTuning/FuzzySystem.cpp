/**
 * @file FuzzySystem.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "FuzzySystem.h"

#include <numeric>

namespace autopas::fuzzy_logic {

void FuzzySystem::addRule(const FuzzyRule &rule) { _rules.push_back(rule); };

std::shared_ptr<FuzzySet> FuzzySystem::applyRules(const std::map<std::string, double> &data) const {
  return std::transform_reduce(_rules.begin(), _rules.end(), _rules[0].apply(data), operator||,
                               [&data](const FuzzyRule &rule) { return rule.apply(data); });
}

double FuzzySystem::predict(const std::map<std::string, double> &data, size_t numSamples) const {
  auto unionSet = applyRules(data);
  return unionSet->centroid(numSamples);
}

}  // namespace autopas::fuzzy_logic