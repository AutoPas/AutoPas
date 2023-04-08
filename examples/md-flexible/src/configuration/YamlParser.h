/**
 * @file YamlParser.h
 * @author N. Fottner
 * @date 15.07.2019
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include "MDFlexConfig.h"
#include "autopas/utils/NumberSet.h"
#include "objects/Object.h"

/**
 * Parser for input through YAML files.
 */
namespace MDFlexParser::YamlParser {

/**
 * Custom Exception for the Yaml parser that is thrown if an error occured during parsing.
 */
class YamlParserException : public std::exception {
 private:
  const char *message;

 public:
  YamlParserException(const char *msg) : message(msg) {}
  const char *what() const noexcept override { return message; }
};

/**
 * Parses the Input for the simulation from the Yaml File specified in the configuration
 * @param config configuration where the input is stored.
 * @return false if any errors occurred during parsing.
 */
bool parseYamlFile(MDFlexConfig &config);
}  // namespace MDFlexParser::YamlParser
