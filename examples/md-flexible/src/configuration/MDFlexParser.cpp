/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#include "MDFlexParser.h"

MDFlexParser::exitCodes MDFlexParser::parseInput(int argc, char **argv, MDFlexConfig &config) {
  // we need to copy argv because the call to getOpt in _cliParser.inputFilesPresent reorders it...
  auto **argvCopy = new char *[argc + 1];
  for (int i = 0; i < argc; i++) {
    auto len = std::string(argv[i]).length() + 1;
    argvCopy[i] = new char[len];
    strcpy(argvCopy[i], argv[i]);
  }
  argvCopy[argc] = nullptr;

  // local helper function to clean up manually allocated arrays
  auto cleanUpArgvCopy = [&]() {
    for (int i = 0; i < argc; i++) {
      delete[] argvCopy[i];
    }
    delete[] argvCopy;
  };

  CLIParser::inputFilesPresent(argc, argv, config);

  MDFlexParser::exitCodes exitCode{MDFlexParser::exitCodes::success};
  if (not config.yamlFilename.value.empty()) {
    // parseYamlFile might throw. Make sure to clean up either way.
    try {
      if (not YamlParser::parseYamlFile(config)) {
        exitCode = MDFlexParser::exitCodes::parsingError;
      }
    } catch (const std::runtime_error &e) {
      cleanUpArgvCopy();
      throw e;
    }
  }

  if (exitCode != MDFlexParser::exitCodes::parsingError) {
    exitCode = CLIParser::parseInput(argc, argvCopy, config);
  }

  cleanUpArgvCopy();

  return exitCode;
}
