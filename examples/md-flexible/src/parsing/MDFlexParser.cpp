/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  // we need to copy argv because the call to getOpt in _cliParser.inputFilesPresent reorders it...
  auto *argvCopy = new char *[argc + 1];
  for (int i = 0; i < argc; i++) {
    auto len = std::string(argv[i]).length() + 1;
    argvCopy[i] = new char[len];
    strncpy(argvCopy[i], argv[i], len);
  }
  argvCopy[argc] = nullptr;

  CLIParser::inputFilesPresent(argc, argv, config);

  if (not config.yamlFilename.value.empty()) {
    if (not YamlParser::parseYamlFile(config)) {
      return false;
    }
  }

  auto parseSuccess = CLIParser::parseInput(argc, argvCopy, config);

  for (int i = 0; i < argc; i++) {
    delete[] argvCopy[i];
  }
  delete[] argvCopy;

  return parseSuccess;
}
