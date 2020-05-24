/**
 * @file CLIParser.h
 * @author F. Gratl
 * @date 18.10.2019
 */

#pragma once

#include <getopt.h>

#include <iostream>

#include "MDFlexConfig.h"

/**
 * Parser for input from the command line.
 */
namespace CLIParser {
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
 * @return false if any errors occurred during parsing.
 */
bool parseInput(int argc, char **argv, MDFlexConfig &config);

/**
 * Creates a file "_md-flexible" containing the completions definitions for zsh.
 * @param cliOptions vector of options to include in the file.
 */
void createZSHCompletionFile(const std::vector<MDFlexConfig::MDFlexOptionInterface> &cliOptions);

};  // namespace CLIParser
