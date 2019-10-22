/**
 * @file YamlParser.h
 * @author N. Fottner
 * @date 15/7/19
 */
#pragma once

#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include "MDFlexConfig.h"
#include "Objects.h"
#include "Thermostat.h"
#include "autopas/autopasIncludes.h"
#include "autopas/utils/NumberSet.h"

/**
 * Parser for input through YAML files.
 */
class YamlParser {
 public:
  /**
   * Parses the Input for the simulation from the Yaml File specified in the configuration
   * @param config configuration where the input is stored.
   * @return false if any errors occurred during parsing.
   * @FIXME: at the moment false is never returned
   */
  bool parseYamlFile(MDFlexConfig &config);

 private:

  /**
   * Check if the passed file exists.
   * @param filename
   * @return true iff the file is found.
   */
  static bool checkFileExists(std::string filename);
};
