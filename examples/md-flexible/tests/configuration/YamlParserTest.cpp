/**
 * @file YamlParserTest.cpp
 * @author N. Fottner
 * @date 02/08/19
 */
#include "YamlParserTest.h"

#include "src/configuration/MDFlexConfig.h"
#include "src/configuration/YamlParser.h"

// @todo: These tests are more related to the MDFlexConfig as to the yaml parser
// See if the tests are required and rename the testfile or move them to the yaml parser

TEST_F(YamlParserTest, calcAutoPasBox) {
  std::vector<std::string> arguments =
    { "md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char* argv[3] = { arguments[0].data(), arguments[1].data(), arguments[2].data() };

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  EXPECT_EQ(configuration.boxMin.value, expectedBoxMin);
  std::array<double, 3> expectedBoxMax = {23, 10, 15.75};
  EXPECT_EQ(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(YamlParserTest, addType) {
  std::vector<std::string> arguments =
    { "md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char* argv[3] = { arguments[0].data(), arguments[1].data(), arguments[2].data() };

  MDFlexConfig configuration(3, argv);

  configuration.addParticleType(0, 1.0, 1.0, 1.0);
  EXPECT_NO_THROW(configuration.addParticleType(0, 1.0, 1.0, 1.0));
  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.5, 1.0, 1.0));
  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.5, 1.1, 1.0));
  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.1, 1.1, 1.1));
  EXPECT_NO_THROW(configuration.addParticleType(1, 2.0, 2.0, 2.0));
  EXPECT_EQ(configuration.massMap.value.at(0), 1.0);
  EXPECT_EQ(configuration.massMap.value.at(1), 2.0);
  EXPECT_EQ(configuration.epsilonMap.value.at(1), 2.0);
}

TEST_F(YamlParserTest, wrongTypeParsingInput) {
  std::vector<std::string> arguments =
    { "md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "incorrectParsingFile.yaml"};

  char* argv[3] = { arguments[0].data(), arguments[1].data(), arguments[2].data() };

  ASSERT_ANY_THROW(MDFlexConfig(3, argv));
}

TEST_F(YamlParserTest, multipleSameObjectParsing) {
  std::vector<std::string> arguments =
    { "md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "multipleSimilarObjects.yaml"};

  char* argv[3] = { arguments[0].data(), arguments[1].data(), arguments[2].data() };

  MDFlexConfig configuration(3, argv);

  ASSERT_EQ(configuration.cubeGridObjects.size(), 2);
  ASSERT_EQ(configuration.cubeGridObjects.at(0).getTypeId(), 0);
  ASSERT_EQ(configuration.cubeGridObjects.at(0).getParticleSpacing(), 0.5);
  ASSERT_EQ(configuration.cubeGridObjects.at(1).getTypeId(), 1);
}
