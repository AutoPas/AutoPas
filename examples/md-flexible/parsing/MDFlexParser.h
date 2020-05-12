/**
 * @file MDFlexParser.cpp
 * @author F. Gratl
 * @date 18.10.2019
 */

#pragma once

#include "CLIParser.h"
#include "MDFlexConfig.h"
#include "YamlParser.h"

namespace MDFlexParser {
bool parseInput(int argc, char **argv, MDFlexConfig &config);
};
