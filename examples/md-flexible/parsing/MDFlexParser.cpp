/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 10/18/19
 */

#include "MDFlexParser.h"

bool MDFlexParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  // we need to copy argv because the call to getOpt in _cliParser.yamlFilePresent reorders it...
  auto argvCopy = new char *[argc + 1];
  for (int i = 0; i < argc; i++) {
    auto len = strlen(argv[i]) + 1;
    argvCopy[i] = new char[len];
    strcpy(argvCopy[i], argv[i]);
  }
  argvCopy[argc] = nullptr;

  if (_cliParser.yamlFilePresent(argc, argv, config)) {
    if (not _yamlParser.parseYamlFile(config)) {
      return false;
    }
  }

  auto parseSuccess = _cliParser.parseInput(argc, argvCopy, config);

  for (int i = 0; i < argc; i++) {
    delete[] argvCopy[i];
  }
  delete[] argvCopy;

  return parseSuccess;
}
