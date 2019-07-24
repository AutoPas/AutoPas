
#include "YamlParserTest.h"


TEST_F(YamlParserTest,calcAutopasBox){
std::array<double, 3> compBoxMin = {-5,-15,-20};
std::array<double, 3> compBoxMax = {19,7.5,9};
EXPECT_EQ(parser.getBoxMin(),compBoxMin);
EXPECT_EQ(parser.getBoxMax(),compBoxMax);
}
