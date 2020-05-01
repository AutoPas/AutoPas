/**
 * @file FullSearchMPIMPITest.cpp
 * @author W. Thieme
 * @date 05/01/20
 */

#include "FullSearchMPITest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(FullSearchMPITest, testSearchSpaceEmpty) {
  autopas::FullSearchMPI fullSearchMPI({});
  EXPECT_TRUE(fullSearchMPI.searchSpaceIsEmpty());
  EXPECT_FALSE(fullSearchMPI.searchSpaceIsTrivial());
  EXPECT_THAT(fullSearchMPI.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(FullSearchMPITest, testSearchSpaceOneOption) {
  autopas::FullSearchMPI fullSearchMPI(
          {autopas::Configuration(autopas::ContainerOption::directSum, 1., autopas::TraversalOption::directSumTraversal,
                                  autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)});
  EXPECT_FALSE(fullSearchMPI.searchSpaceIsEmpty());
  EXPECT_TRUE(fullSearchMPI.searchSpaceIsTrivial());
  EXPECT_THAT(fullSearchMPI.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(FullSearchMPITest, testSearchSpaceMoreOptions) {
  autopas::FullSearchMPI fullSearchMPI({autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08},
                                 {autopas::DataLayoutOption::soa},
                                 {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  EXPECT_FALSE(fullSearchMPI.searchSpaceIsEmpty());
  EXPECT_FALSE(fullSearchMPI.searchSpaceIsTrivial());
  EXPECT_THAT(fullSearchMPI.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(FullSearchMPITest, testRemoveN3OptionRemoveAll) {
  autopas::FullSearchMPI fullSearchMPI(
          {autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
          {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

  EXPECT_THROW(fullSearchMPI.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(FullSearchMPITest, testRemoveN3OptionRemoveSome) {
  autopas::FullSearchMPI fullSearchMPI({autopas::ContainerOption::linkedCells}, {1.},
                                 {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                 {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                 {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

  EXPECT_NO_THROW(fullSearchMPI.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(fullSearchMPI.searchSpaceIsEmpty());
  EXPECT_FALSE(fullSearchMPI.searchSpaceIsTrivial());
}

TEST_F(FullSearchMPITest, testTune) {
  autopas::FullSearchMPI fullSearchMPI(
          {autopas::ContainerOption::linkedCells}, {1.},
          {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
          {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearchMPI.getCurrentConfiguration());
  fullSearchMPI.addEvidence(10);

  fullSearchMPI.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearchMPI.getCurrentConfiguration());
  fullSearchMPI.addEvidence(1);

  fullSearchMPI.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearchMPI.getCurrentConfiguration());
  fullSearchMPI.addEvidence(20);

  fullSearchMPI.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearchMPI.getCurrentConfiguration());
}