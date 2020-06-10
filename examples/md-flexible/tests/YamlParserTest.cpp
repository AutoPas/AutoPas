/**
 * @file YamlParserTest.cpp
 * @author N. Fottner
 * @date 02/08/19
 */
#include "YamlParserTest.h"

#include "src/parsing/MDFlexConfig.h"
#include "src/parsing/YamlParser.h"

TEST_F(YamlParserTest, calcAutoPasBox) {
  // tests the AutoBox calculation after Object initialization with multipleObjectsWithMultipleTypeTest
  MDFlexConfig config;
  config.yamlFilename.value = std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
  MDFlexParser::YamlParser::parseYamlFile(config);
  config.calcSimulationBox();
  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  EXPECT_EQ(config.boxMin.value, expectedBoxMin);
  std::array<double, 3> expectedBoxMax = {23, 10, 15.75};
  EXPECT_EQ(config.boxMax.value, expectedBoxMax);
}

TEST_F(YamlParserTest, addType) {
  // tests if exception are thrown if types arent well initialized
  MDFlexConfig config;
  config.addParticleType(0, 1.0, 1.0, 1.0);
  EXPECT_NO_THROW(config.addParticleType(0, 1.0, 1.0, 1.0));
  EXPECT_ANY_THROW(config.addParticleType(0, 1.5, 1.0, 1.0));
  EXPECT_ANY_THROW(config.addParticleType(0, 1.5, 1.1, 1.0));
  EXPECT_ANY_THROW(config.addParticleType(0, 1.1, 1.1, 1.1));
  EXPECT_NO_THROW(config.addParticleType(1, 2.0, 1.0, 1.0));
  EXPECT_EQ(config.massMap.value.at(0), 1.0);
  EXPECT_EQ(config.massMap.value.at(1), 1.0);
  EXPECT_EQ(config.epsilonMap.value.at(1), 2.0);
}

TEST_F(YamlParserTest, wrongTypeParsingInput) {
  MDFlexConfig config;
  config.yamlFilename.value = std::string(YAMLDIRECTORY) + "incorrectParsingFile.yaml";
  ASSERT_ANY_THROW(MDFlexParser::YamlParser::parseYamlFile(config));
}

TEST_F(YamlParserTest, multipleSameObjectParsing) {
  MDFlexConfig config;
  config.yamlFilename.value = std::string(YAMLDIRECTORY) + "multipleSimilarObjects.yaml";
  MDFlexParser::YamlParser::parseYamlFile(config);
  ASSERT_EQ(config.cubeGridObjects.size(), 2);
  ASSERT_EQ(config.cubeGridObjects.at(0).getTypeId(), 0);
  ASSERT_EQ(config.cubeGridObjects.at(0).getParticleSpacing(), 0.5);
  ASSERT_EQ(config.cubeGridObjects.at(1).getTypeId(), 1);
}