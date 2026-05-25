/**
 * @file OutputMapper.h
 * @author Manuel Lerchner
 * @date 29.05.24
 */

#pragma once

#include <map>
#include <memory>
#include <string>

#include "autopas/tuning/tuningStrategy/ruleBasedTuning/RuleBasedProgramTree.h"

namespace autopas::FuzzyLogic {

/**
 * Used to store the mapping information to transform the result of the fuzzy controller into a configuration suitable
 * for AutoPas.
 */
class OutputMapper {
 public:
  /**
   * Constructor of the OutputMapper.
   * @param outputDomain The domain of the output variable.
   * @param mappings A vector of pairs of the form {value, configurationPatterns}.
   */
  OutputMapper(std::string outputDomain, std::vector<std::pair<double, std::vector<ConfigurationPattern>>> mappings);

  /**
   * Returns a string representation of the OutputMapping
   */
  explicit operator std::string() const;

  /**
   * Returns the closest configuration patterns to the given value.
   * @param value The value to search for.
   * @return The closest configuration patterns to the given value.
   */
  [[nodiscard]] std::vector<ConfigurationPattern> getClosestConfigurationPatterns(double value);

  /**
   * Getter for the output domain.
   * @return The output domain.
   */
  [[nodiscard]] const std::string &getOutputDomain();

  /**
   * Getter for the mappings.
   * @return The mappings.
   */
  [[nodiscard]] const std::vector<std::pair<double, std::vector<ConfigurationPattern>>> &getMappings() const {
    return _mappings;
  }

 private:
  /**
   * The name of the domain which should be mapped.
   */
  std::string _outputDomain;

  /**
   * The mappings from the fuzzy output to the configuration patterns.
   */
  std::vector<std::pair<double, std::vector<ConfigurationPattern>>> _mappings;
};

}  // namespace autopas::FuzzyLogic
