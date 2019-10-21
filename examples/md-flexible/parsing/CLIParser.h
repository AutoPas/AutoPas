/**
 * @file CLIParser.h
 * @author F. Gratl
 * @date 10/18/19
 */

#pragma once

#include <getopt.h>
#include "MDFlexConfig.h"
#include <iostream>

/**
 * Parser for input from the command line.
 */
class CLIParser {
 public:
  /**
   * Checks if a yaml file is specified in the given command line arguments and saves it's path to the given
   * configuration.
   * @param argc number of command line arguments.
   * @param argv command line argument array.
   * @param config configuration where the input is stored.
   * @return true if a yaml file argument was found.
   */
  bool yamlFilePresent(int argc, char **argv, MDFlexConfig &config);

  /**
   * Parses the Input for the simulation from the command line.
   * @param argc number of command line arguments.
   * @param argv command line argument array.
   * @param config configuration where the input is stored.
   * @return false if any errors occurred during parsing.
   */
  bool parseInput(int argc, char **argv, MDFlexConfig &config);
};

