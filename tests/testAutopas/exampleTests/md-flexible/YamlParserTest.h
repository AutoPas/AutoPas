#pragma once
#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include "Objects/Objects.h"
#include "../../../../examples/md-flexible/parsing/MDFlexConfig.h"
#include "../../../../examples/md-flexible/parsing/YamlParser.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/commonTypedefs.h"

class YamlParserTest : public AutoPasTestBase {
 public:
  YamlParserTest() : AutoPasTestBase() {}
};
