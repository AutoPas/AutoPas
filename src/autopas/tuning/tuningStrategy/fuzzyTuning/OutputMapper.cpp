/**
 * @file OutputMapper.cpp
 * @author Manuel Lerchner
 * @date 29.05.24
 */

#include "OutputMapper.h"

autopas::fuzzy_logic::OutputMapper::OutputMapper(
    std::string outputDomain, std::vector<std::pair<double, std::vector<ConfigurationPattern>>> mappings)
    : _outputDomain(std::move(outputDomain)), _mappings(std::move(mappings)) {}

const std::string &autopas::fuzzy_logic::OutputMapper::getOutputDomain() { return _outputDomain; }

autopas::fuzzy_logic::OutputMapper::operator std::string() const {
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
