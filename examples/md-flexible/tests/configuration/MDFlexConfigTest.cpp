/**
 * @file MDFlexConfigTest.cpp
 * @author F. Gratl
 * @date 04/06/2020
 */

#include "MDFlexConfigTest.h"

#include <gmock/gmock-matchers.h>

using ::testing::ElementsAreArray;

TEST_F(MDFlexConfigTest, GridBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-file", std::string(YAMLDIRECTORY) + "cubeGrid.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {0, 0, 0};
  std::array<double, 3> expectedBoxMax = {9.75, 9.75, 9.75};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, SphereBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "sphere.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  std::array<double, 3> expectedBoxMax = {20.75, 5.75, 15.75};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, GaussBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "cubeGauss.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {0.0, 0.0, 0.0};
  std::array<double, 3> expectedBoxMax = {23.0, 8.0, 13.0};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, UniformBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "cubeUniform.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {0, 0, -15};
  std::array<double, 3> expectedBoxMax = {10, 10, 1};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, ClosestPackedBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "cubeClosestPacked.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {-0.5, -0.5, 0.};
  std::array<double, 3> expectedBoxMax = {7., 7.5, 27.5};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, calcAutoPasBox) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  EXPECT_EQ(configuration.boxMin.value, expectedBoxMin);
  std::array<double, 3> expectedBoxMax = {23, 10, 27.5};
  EXPECT_EQ(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, addType) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

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

TEST_F(MDFlexConfigTest, wrongTypeParsingInput) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "incorrectParsingFile.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  // If an invalid YAML-file is used, exceptions are catched by YamlParser, MDFlexConfig will then exit with
  // EXIT_FAILURE and write "Error when parsing configuration file." to cerr. YAML-file for this test is in
  // examples/md-flexible/tests/yamlTestFiles/incorrectParsingFile.yaml
  ASSERT_EXIT(MDFlexConfig(3, argv), testing::ExitedWithCode(1), "Error when parsing configuration file.");
}

TEST_F(MDFlexConfigTest, multipleSameObjectParsing) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "multipleSimilarObjects.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  ASSERT_EQ(configuration.cubeGridObjects.size(), 2);
  ASSERT_EQ(configuration.cubeGridObjects.at(0).getTypeId(), 0);
  ASSERT_EQ(configuration.cubeGridObjects.at(0).getParticleSpacing(), 0.5);
  ASSERT_EQ(configuration.cubeGridObjects.at(1).getTypeId(), 1);
}
