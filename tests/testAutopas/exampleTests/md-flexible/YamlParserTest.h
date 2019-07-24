#pragma  once
#include <gtest/gtest.h>
#include <math.h>
#include <vector>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"
#include "../../../../examples/md-flexible/Objects.h"
#include "../../../../examples/md-flexible/YamlParser.h"
#include "autopas/utils/ArrayUtils.h"


class YamlParserTest  :public AutoPasTestBase {
public:
    YamlParserTest() :AutoPasTestBase() ,parser{YamlParser()} {parser.parseInput(filename);}

protected:
    YamlParser parser;
    std::string filename="testParsing.yaml";
};
