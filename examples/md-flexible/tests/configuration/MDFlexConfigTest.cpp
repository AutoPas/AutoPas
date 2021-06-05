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

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {0, 0, 0};
  std::array<double, 3> expectedBoxMax = {9.75, 9.75, 9.75};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, SphereBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "sphere.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {-10.75, -25.75, -15.75};
  std::array<double, 3> expectedBoxMax = {20.75, 5.75, 15.75};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, GaussBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "cubeGauss.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {0.0, 0.0, 0.0};
  std::array<double, 3> expectedBoxMax = {23.0, 8.0, 13.0};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}

TEST_F(MDFlexConfigTest, UniformBoxMinMax) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "cubeUniform.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> expectedBoxMin = {0, 0, -15};
  std::array<double, 3> expectedBoxMax = {10, 10, 1};

  EXPECT_THAT(configuration.boxMin.value, expectedBoxMin);
  EXPECT_THAT(configuration.boxMax.value, expectedBoxMax);
}
