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
	std::string arguments = "md-flexible --yaml-filename "
		+ std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";

  MDFlexConfig configuration(3, reinterpret_cast<char**>(arguments.data()));
  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  EXPECT_EQ(configuration.boxMin.value, expectedBoxMin);
  std::array<double, 3> expectedBoxMax = {23, 10, 15.75};
  EXPECT_EQ(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(YamlParserTest, addType) {
	std::string arguments = "md-flexible --yaml-filename "
		+ std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml";
  MDFlexConfig configuration(3, reinterpret_cast<char**>(arguments.data()));

  configuration.addParticleType(0, 1.0, 1.0, 1.0);
  EXPECT_NO_THROW(configuration.addParticleType(0, 1.0, 1.0, 1.0));
  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.5, 1.0, 1.0));
  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.5, 1.1, 1.0));
  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.1, 1.1, 1.1));
  EXPECT_NO_THROW(configuration.addParticleType(1, 2.0, 1.0, 1.0));
  EXPECT_EQ(configuration.massMap.value.at(0), 1.0);
  EXPECT_EQ(configuration.massMap.value.at(1), 1.0);
  EXPECT_EQ(configuration.epsilonMap.value.at(1), 2.0);
}

TEST_F(YamlParserTest, wrongTypeParsingInput) {
	std::string arguments = "md-flexible --yaml-filename "
		+ std::string(YAMLDIRECTORY) + "incorrectParsingFile.yaml";

  ASSERT_ANY_THROW(MDFlexConfig(3, reinterpret_cast<char**>(arguments.data())));
}

TEST_F(YamlParserTest, multipleSameObjectParsing) {
	std::string arguments = "md-flexible --yaml-filename "
		+ std::string(YAMLDIRECTORY) + "multipleSimilarObjects.yaml";
  MDFlexConfig configuration(3, reinterpret_cast<char**>(arguments.data()));

  ASSERT_EQ(configuration.cubeGridObjects.size(), 2);
  ASSERT_EQ(configuration.cubeGridObjects.at(0).getTypeId(), 0);
  ASSERT_EQ(configuration.cubeGridObjects.at(0).getParticleSpacing(), 0.5);
  ASSERT_EQ(configuration.cubeGridObjects.at(1).getTypeId(), 1);
}
