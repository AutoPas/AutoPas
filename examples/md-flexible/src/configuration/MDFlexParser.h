/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#pragma once

#include "CLIParser.h"
#include "MDFlexConfig.h"
#include "ParserExitCodes.h"
#include "YamlParser.h"

/**
 * General parser.
 */
namespace MDFlexParser {

/**
 * Parse the given command line options and if necessary also the yaml file within that.
 * Parsed values are directly stored to the passed config object.
 * @param argc
 * @param argv
 * @param config
 * @return Indicator of success. See MDFlexParser::exitCodes for possible values.
 */
exitCodes parseInput(int argc, char **argv, MDFlexConfig &config);

}
