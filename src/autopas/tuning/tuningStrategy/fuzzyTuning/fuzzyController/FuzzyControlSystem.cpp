/**
 * @file FuzzyControlSystem.cpp
 * @author Manuel Lerchner
 * @date 18.04.24
 */

#include "FuzzyControlSystem.h"

#include <numeric>

#include "autopas/utils/ExceptionHandler.h"

namespace autopas::FuzzyLogic {

FuzzyControlSystem::FuzzyControlSystem(std::shared_ptr<FuzzyControlSettings> settings) : _settings(settings) {}

void FuzzyControlSystem::addRule(const FuzzyRule &rule) {
  auto dimensions = rule.getConsequent()->getCrispSet()->getDimensions();
  if (dimensions.size() != 1) {
    autopas::utils::ExceptionHandler::exception("The consequent of a rule must be a one-dimensional fuzzy set");
  }
  auto [dimensionName, _] = *dimensions.begin();

  if (not _outputDomain.has_value()) {
    _outputDomain = dimensionName;
  } else if (_outputDomain.value() != dimensionName) {
    autopas::utils::ExceptionHandler::exception("All consequents of the FuzzyControlSystem must have the same domain");
  }

  _rules.push_back(rule);
};

std::shared_ptr<FuzzySet> FuzzyControlSystem::applyRules(const FuzzySet::Data &data) const {
  // Applies all the rules of the system to get individual fuzzy sets corresponding to the cut-consequents.
  // Then combines them using the Fuzzy-OR operator to get the union set.
  return std::transform_reduce(_rules.begin() + 1, _rules.end(), _rules[0].apply(data), operator||,
                               [&data](const FuzzyRule &rule) { return rule.apply(data); });
}

double FuzzyControlSystem::predict(const FuzzySet::Data &data, size_t numSamples) const {
  auto unionSet = applyRules(data);

  if (_settings->count("defuzzificationMethod") == 0) {
    autopas::utils::ExceptionHandler::exception("No defuzzification method specified in the settings");
  }

  if (_settings->count("numSamples") != 0) {
    numSamples = std::stoul(_settings->at("numSamples"));
  }

  auto defuzzificationMethod = DefuzzificationMethodOption::parseOptionExact(_settings->at("defuzzificationMethod"));

  AutoPasLog(DEBUG, "Defuzzifying with method: {} and numSamples: {}", defuzzificationMethod.to_string(), numSamples);

  return unionSet->defuzzify(defuzzificationMethod, numSamples);
}

FuzzyControlSystem::operator std::string() const {
  return std::accumulate(_rules.begin(), _rules.end(),
                         std::string("FuzzyControlSystem: \"" + _outputDomain.value() + "\"\n"),
                         [](const std::string &acc, const FuzzyRule &rule) { return acc + std::string(rule) + "\n"; });
}

}  // namespace autopas::FuzzyLogic