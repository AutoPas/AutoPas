#include "YamlParserTest.h"

TEST_F(YamlParserTest, calcAutopasBox) {
  //tests the AutoBox calculation after Object initialization with multipleObjectsWithMultipleTypeTest
    std::string file = "multipleObjectsWithMultipleTypesTest.yaml";
    parser.setFilename(file);
    parser.parseYamlFile();
    std::array<double, 3> compBoxMin = {-5, -15, -20};
  std::array<double, 3> compBoxMax = {19, 7.5, 9};
  EXPECT_EQ(parser.getBoxMin(), compBoxMin);
  EXPECT_EQ(parser.getBoxMax(), compBoxMax);
}

TEST_F(YamlParserTest,addType){
    //tests if exception are thrown if types arent well initialized
    std::map<unsigned long, double> epsilonMap;
    std::map<unsigned long, double> sigmaMap ;
    std::map<unsigned long, double> massMap;
    parser.addType(0,1.0,1.0,1.0);
    EXPECT_NO_THROW(parser.addType(0,1.0,1.0,1.0));
    EXPECT_ANY_THROW(parser.addType(0,1.5,1.0,1.0));
    EXPECT_ANY_THROW(parser.addType(0,1.5,1.1,1.0));
    EXPECT_ANY_THROW(parser.addType(0,1.1,1.1,1.1));
    EXPECT_NO_THROW(parser.addType(1,2.0,1.0,1.0));
    EXPECT_EQ(parser.getMassMap().at(0),1.0);
    EXPECT_EQ(parser.getMassMap().at(1),1.0);
    EXPECT_EQ(parser.getEpsilonMap().at(1),2.0);
}

TEST_F(YamlParserTest,wrongTypeParsingInput){
    std::string file = "incorrectParsingFile.yaml";
    this->parser.setFilename(file);
    ASSERT_ANY_THROW(this->parser.parseYamlFile());
}

