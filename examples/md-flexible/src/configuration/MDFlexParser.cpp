/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#include "MDFlexParser.h"

MDFlexParser::exitCodes MDFlexParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  CLIParser::inputFilesPresent(argc, argv, config);

  if (not config.yamlFilename.value.empty()) {
    if (not YamlParser::parseYamlFile(config)) {
      return MDFlexParser::exitCodes::parsingError;
    }
  }

  auto cliParseExitCode = CLIParser::parseInput(argc, argv, config);

  return cliParseExitCode;
}
