/**
 * @file YamlParserTest.h
 * @author N. Fottner
 * @date 02/08/19
 */
#pragma once
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "parsing/MDFlexConfig.h"
#include "parsing/YamlParser.h"
#include "AutoPasTestBase.h"
#include "Objects/Objects.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/commonTypedefs.h"

class YamlParserTest : public AutoPasTestBase {
 public:
  YamlParserTest() : AutoPasTestBase() {}
};
