/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#include "MDFlexConfig.h"
#include "YamlParser.h"
#include "CLIParser.h"

class MDFlexParser {

 public:

  bool parseInput(int argc, char **argv);

  const MDFlexConfig &getConfig() const;

 private:

  CLIParser _cliParser;
  YamlParser _yamlParser;
  MDFlexConfig _config;
};
