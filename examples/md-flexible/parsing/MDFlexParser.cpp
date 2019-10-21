/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#include "MDFlexParser.h"

const MDFlexConfig &MDFlexParser::getConfig() const { return _config; }

bool MDFlexParser::parseInput(int argc, char **argv) {
  if (_cliParser.yamlFilePresent(argc, argv, _config)) {
    if (not _yamlParser.parseYamlFile(_config)){
      return false;
    }
  }

  return _cliParser.parseInput(argc, argv, _config);
}
