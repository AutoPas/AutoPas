/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#pragma once

#include "CLIParser.h"
#include "MDFlexConfig.h"
#include "YamlParser.h"

class MDFlexParser {
 public:
  bool parseInput(int argc, char **argv, MDFlexConfig &config);

 private:
  CLIParser _cliParser;
  YamlParser _yamlParser;
};
