/**
 * @file OutputMapper.cpp
 * @author Manuel Lerchner
 * @date 29.05.24
 */

#include "OutputMapper.h"

using namespace autopas::fuzzy_logic;

OutputMapper::OutputMapper(std::string outputDomain,
                           std::vector<std::pair<double, std::vector<ConfigurationPattern>>> mappings)
    : _outputDomain(std::move(outputDomain)), _mappings(std::move(mappings)) {}

std::vector<autopas::ConfigurationPattern> OutputMapper::getClosestConfigurationPatterns(double value) {
  auto it = std::min_element(_mappings.begin(), _mappings.end(),
                             [value](const std::pair<double, std::vector<ConfigurationPattern>> &a,
                                     const std::pair<double, std::vector<ConfigurationPattern>> &b) {
                               return std::abs(a.first - value) < std::abs(b.first - value);
                             });

  return it->second;
}

const std::string &OutputMapper::getOutputDomain() { return _outputDomain; }

OutputMapper::operator std::string() const {
  std::string result = "OutputMapper: " + _outputDomain + "\n";
  for (const auto &mapping : _mappings) {
    result += "  " + std::to_string(mapping.first) + " -> ";
    for (const auto &pattern : mapping.second) {
      result += "[" + pattern.toString() + "] ";
    }
    result += "\n";
  }
  return result;
}
