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

TEST_F(MDFlexConfigTest, calcAutoPasBox) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};

  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  EXPECT_EQ(configuration.boxMin.value, expectedBoxMin);
  std::array<double, 3> expectedBoxMax = {23, 10, 15.75};
  EXPECT_EQ(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, addType) {
  GTEST_SKIP_("This test needs adapting to multi-site molecules");
//
//  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
//                                        std::string(YAMLDIRECTORY) + "multipleObjectsWithMultipleTypesTest.yaml"};
//
//  char *argv[3] = {&arguments[0][0], &arguments[1][0], &arguments[2][0]};
//
//  MDFlexConfig configuration(3, argv);
//
//  configuration.addParticleType(0, 1.0, 1.0, 1.0);
//  EXPECT_NO_THROW(configuration.addParticleType(0, 1.0, 1.0, 1.0));
//  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.5, 1.0, 1.0));
//  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.5, 1.1, 1.0));
//  EXPECT_ANY_THROW(configuration.addParticleType(0, 1.1, 1.1, 1.1));
//  EXPECT_NO_THROW(configuration.addParticleType(1, 2.0, 2.0, 2.0));
//  EXPECT_EQ(configuration.massMap.value.at(0), 1.0);
//  EXPECT_EQ(configuration.massMap.value.at(1), 2.0);
//  EXPECT_EQ(configuration.epsilonMap.value.at(1), 2.0);
}

TEST_F(MDFlexConfigTest, multipleSimilarObjectParsing) {
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
                                        std::string(YAMLDIRECTORY) + "siteParsingTestCorrect.yaml"};

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
 * Test for the correct parsing of site types from a .yaml file.
 *
 * Loads two incorrect files and checks that errors are thrown. In the first, a non-consecutive site number is given.
 * In the second, incomplete site information is given.
 */
TEST_F(MDFlexConfigTest, incorrectSiteParsing) {
  // Test for non-consecutive site id.
  std::vector<std::string> argumentsTest1 = {"md-flexible", "--yaml-filename",
                                             std::string(YAMLDIRECTORY) + "siteParsingTestIncorrect1.yaml"};

  char *argv1[3] = {&argumentsTest1[0][0], &argumentsTest1[1][0], &argumentsTest1[2][0]};

  // Test error thrown.
  EXPECT_ANY_THROW(MDFlexConfig configuration(3, argv1););

  // Test for incomplete site information.
  std::vector<std::string> argumentsTest2 = {"md-flexible", "--yaml-filename",
                                            std::string(YAMLDIRECTORY) + "siteParsingTestIncorrect2.yaml"};

  char *argv2[3] = {&argumentsTest2[0][0], &argumentsTest2[1][0], &argumentsTest2[2][0]};

  // Test error thrown.
  EXPECT_ANY_THROW(MDFlexConfig configuration(3, argv2););
}

/**
 * Test for the correct parsing of molecule type from a .yaml file.
 */
TEST_F(MDFlexConfigTest, correctMolParsing) {
  // Configure md-flexible using correct yaml file and expect no throw in doing so.
  std::vector<std::string> argumentsTest = {"md-flexible", "--yaml-filename",
                                            std::string(YAMLDIRECTORY) + "molParsingTestCorrect.yaml"};

  char *argv[3] = {&argumentsTest[0][0], &argumentsTest[1][0], &argumentsTest[2][0]};

  MDFlexConfig configuration(3, argv);

  // Expected molecular information
  const std::vector<int> expectedSiteIds0 = {0};
  const std::vector<std::array<double, 3>> expectedSitePositions0 = {{0., 0., 0.}};
  const std::array<double, 3> expectedMoI0 = {1., 1., 1.};
  const std::vector<int> expectedSiteIds1 = {0, 0, 1};
  const std::vector<std::array<double, 3>> expectedSitePositions1 = {{0., -0.5, 0.}, {1., 0.1, 0.01}, {-0.2, -0.3, 0.4}};
  const std::array<double, 3> expectedMoI1 = {1., 0.5, 0.25};

  const auto siteIds0 = configuration.molToSiteIdMap.value.at(0);
  for (int i = 0; i < expectedSiteIds0.size(), i++) {
    EXPECT_EQ(expectedSiteIds0[i], siteIds0) << "Molecule 0's site IDs do not match expected value. Expected:
  }

}