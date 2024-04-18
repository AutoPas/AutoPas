/**
 * @file FuzzySystem.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#pragma once

#include "FuzzySystem.h"

namespace autopas::fuzzy_logic {

void FuzzySystem::addRule(const FuzzyRule &rule) { _rules.push_back(rule); };

FuzzySet FuzzySystem::applyRules(const std::map<std::string, double> &data) const {
  std::vector<FuzzySet> consequents;
  for (const auto &rule : _rules) {
    consequents.push_back(rule.apply(data));
  }

  FuzzySet result = consequents[0];
  for (size_t i = 1; i < consequents.size(); ++i) {
    result = result || consequents[i];
  }

  return result;
}

double FuzzySystem::predict(const std::map<std::string, double> &data, size_t numSamples) const {
  auto unionSet = applyRules(data);
  return unionSet.centroid(numSamples);
}

}  // namespace autopas::fuzzy_logic