/**
 * @file CLIParser.h
 * @author F. Gratl
 * @date 18.10.2019
 */

#pragma once

#include <getopt.h>

#include <iostream>

#include "MDFlexConfig.h"
#include "ParserExitCodes.h"

/**
 * Parser for input from the command line.
 */
namespace MDFlexParser::CLIParser {
/**
 * Checks if a checkpoint or yaml file is specified in the given command line arguments.
 * If any files are found, their paths are saves in the respective fields of the given configuration.
 *
 * @param argc number of command line arguments.
 * @param argv command line argument array.
 * @param config configuration where the input is stored.
 * @return
 */
void inputFilesPresent(int argc, char **argv, MDFlexConfig &config);

/**
 * Parses the Input for the simulation from the command line.
 * @param argc number of command line arguments.
 * @param argv command line argument array.
 * @param config configuration where the input is stored.
 * @return Indicator of success. See MDFlexParser::exitCodes for possible values.
 */
MDFlexParser::exitCodes parseInput(int argc, char **argv, MDFlexConfig &config);

/**
 * Prints the help message to the given stream.
 * @param ostream Typically std::out.
 * @param relPathOfExecutable  Typically argv[0].
 * @param relevantOptions vector of options to include in the help message.
 */
void printHelpMessage(std::ostream &ostream, const std::string &relPathOfExecutable,
                      const std::vector<MDFlexConfig::MDFlexOptionInterface> &relevantOptions);

/**
 * Creates a file "_md-flexible" containing the completions definitions for zsh.
 * The function expects possible options to be in the description of an option and in the form of:
 * 'Possible Values: (val1 val2 ... valN)'
 * @param cliOptions vector of options to include in the file.
 */
void createZSHCompletionFile(const std::vector<MDFlexConfig::MDFlexOptionInterface> &cliOptions);

};  // namespace MDFlexParser::CLIParser
