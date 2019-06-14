/**
 * @file MachineSearchTest.cpp
 * @date 7/2/19
 */

#include "MachineSearchTest.h"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include "testingHelpers/commonTypedefs.h"

TEST_F(MachineSearchTest, testSearchSpaceExpectedOptions) {
  autopas::MachineSearch<Particle, FPCell> machineSearch(
      {autopas::ContainerOption::linkedCells, autopas::ContainerOption::verletLists,
       autopas::ContainerOption::verletListsCells}, {1.},
      autopas::allTraversalOptions, autopas::allDataLayoutOptions, autopas::allNewton3Options, "fdeep_model.json");
  EXPECT_FALSE(machineSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(machineSearch.searchSpaceIsTrivial());
  EXPECT_THAT(machineSearch.getAllowedContainerOptions(), ::testing::Contains(autopas::ContainerOption::linkedCells));
}

TEST_F(MachineSearchTest, testTune) {
  autopas::MachineSearch<Particle, FPCell> machineSearch(
      {autopas::ContainerOption::linkedCells},{1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, "fdeep_model.json");

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            machineSearch.getCurrentConfiguration());
  machineSearch.addEvidence(10);

  machineSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            machineSearch.getCurrentConfiguration());
  machineSearch.addEvidence(1);

  machineSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            machineSearch.getCurrentConfiguration());
  machineSearch.addEvidence(20);

  machineSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            machineSearch.getCurrentConfiguration());
}