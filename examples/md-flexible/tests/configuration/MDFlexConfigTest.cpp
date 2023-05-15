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

// todo: add test for parsing site and molecule information

/**
 * Test for the correct parsing of site types from a .yaml file.
 */
TEST_F(MDFlexConfigTest, correctSiteParsing) {
  // Configure md-flexible using correct yaml file and expect no throw in doing so.
  std::vector<std::string> argumentsTest = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "siteParsingTest.yaml"};

  char *argv[3] = {&argumentsTest[0][0], &argumentsTest[1][0], &argumentsTest[2][0]};

  MDFlexConfig configuration(3, argv);

  // Test correct site information obtained.
  EXPECT_EQ(configuration.epsilonMap.value.at(0), 1.);
  EXPECT_EQ(configuration.sigmaMap.value.at(0), 1.);
  EXPECT_EQ(configuration.massMap.value.at(0), 1.);
  EXPECT_EQ(configuration.epsilonMap.value.at(1), 1.2);
  EXPECT_EQ(configuration.sigmaMap.value.at(1), 0.8);
  EXPECT_EQ(configuration.massMap.value.at(1), 0.5);
}

/**
 * Test for the correct parsing of molecule type from a .yaml file.
 */
TEST_F(MDFlexConfigTest, correctMolParsing) {
  // Skip test if not compiled for multi-site molecules.
#if not MD_FLEXIBLE_MODE==MULTISITE
  GTEST_SKIP() << "correctMolParsing: Skipping as multi-site not compiled";
#endif

  // Configure md-flexible using correct yaml file and expect no throw in doing so.
  std::vector<std::string> argumentsTest = {"md-flexible", "--yaml-filename",
                                            std::string(YAMLDIRECTORY) + "molParsingTest.yaml"};

  char *argv[3] = {&argumentsTest[0][0], &argumentsTest[1][0], &argumentsTest[2][0]};

  MDFlexConfig configuration(3, argv);

  // Expected molecular information
  const std::vector<int> expectedSiteIds0 = {0};
  const std::vector<std::array<double, 3>> expectedSitePositions0 = {{0., 0., 0.}};
  const std::array<double, 3> expectedMoI0 = {1., 1., 1.};
  const std::vector<int> expectedSiteIds1 = {0, 0, 1};
  const std::vector<std::array<double, 3>> expectedSitePositions1 = {{0., -0.5, 0.}, {1., 0.1, 0.01}, {-0.2, -0.3, 0.4}};
  const std::array<double, 3> expectedMoI1 = {1., 0.5, 0.25};

  const auto siteIds0 = configuration.molToSiteIdMap.at(0);
  for (int i = 0; i < expectedSiteIds0.size(); i++) {
    EXPECT_EQ(expectedSiteIds0[i], siteIds0[i]);
  }
  const auto sitePositions0 = configuration.molToSitePosMap.at(0);
  for (int i = 0; i < expectedSitePositions0.size(); i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(expectedSitePositions0[i][j], sitePositions0[i][j]);
    }
  }
  const auto momentOfInertia0 = configuration.momentOfInertiaMap.at(0);
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(expectedMoI0[i], momentOfInertia0[i]);
  }

  const auto siteIds1 = configuration.molToSiteIdMap.at(1);
  for (int i = 0; i < expectedSiteIds1.size(); i++) {
    EXPECT_EQ(expectedSiteIds1[i], siteIds1[i]);
  }
  const auto sitePositions1 = configuration.molToSitePosMap.at(1);
  for (int i = 0; i < expectedSitePositions1.size(); i++) {
    for (int j = 0; j < 3; j++) {
      EXPECT_EQ(expectedSitePositions1[i][j], sitePositions1[i][j]);
    }
  }
  const auto momentOfInertia1 = configuration.momentOfInertiaMap.at(1);
  for (int i = 0; i < 3; i++) {
    EXPECT_EQ(expectedMoI1[i], momentOfInertia1[i]);
  }
}