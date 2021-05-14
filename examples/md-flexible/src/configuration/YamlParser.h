/**
 * @file YamlParser.h
 * @author N. Fottner
 * @date 15.07.2019
 */
#pragma once

#include <yaml-cpp/yaml.h>

#include "autopas/utils/NumberSet.h"
#include "MDFlexConfig.h"
#include "objects/Object.h"

/**
 * Parser for input through YAML files.
 */
namespace MDFlexParser::YamlParser {
	/**
 	* Parses the Input for the simulation from the Yaml File specified in the configuration
 	* @param config configuration where the input is stored.
 	* @return false if any errors occurred during parsing.
 	* @note FIXME: at the moment false is never returned and the parser just ungracefully crashes.
 	*/
	bool parseYamlFile(MDFlexConfig &config);
}
