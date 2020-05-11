/**
 * @file FullSearchMPIMPITest.cpp
 * @author W. Thieme
 * @date 05/01/20
 */

#include "FullSearchMPITest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <unistd.h>

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
  EXPECT_THAT(fullSearchMPI.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(FullSearchMPITest, testRemoveN3OptionRemoveAll) {
  autopas::FullSearchMPI fullSearchMPI(
          {autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
          {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

  EXPECT_THROW(fullSearchMPI.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(FullSearchMPITest, testGlobalOptimumAndReset) {
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);

  autopas::FullSearchMPI fullSearchMPI(
          {autopas::Configuration(autopas::ContainerOption::directSum, 1. + (double)worldRank/10.,
                                  autopas::TraversalOption::directSumTraversal,
                                  autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)});
  fullSearchMPI.addEvidence(worldRank);

  fullSearchMPI.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::directSum, 1. + (double) worldRank/10.,
                                   autopas::TraversalOption::directSumTraversal,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());

  // test synchronization
  usleep(worldRank * 1000000);
  fullSearchMPI.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::directSum, 1., autopas::TraversalOption::directSumTraversal,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());

  fullSearchMPI.reset();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::directSum, 1. + (double)worldRank/10.,
                                   autopas::TraversalOption::directSumTraversal,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());
}

TEST_F(FullSearchMPITest, testLocalOptimumAndReset) {
  autopas::FullSearchMPI fullSearchMPI({autopas::Configuration(autopas::ContainerOption::directSum, 1.,
                                                               autopas::TraversalOption::directSumTraversal,
                                                               autopas::DataLayoutOption::soa,
                                                               autopas::Newton3Option::enabled),
                                        autopas::Configuration(autopas::ContainerOption::linkedCells, 1.2,
                                                               autopas::TraversalOption::c08,
                                                               autopas::DataLayoutOption::aos,
                                                               autopas::Newton3Option::enabled),
                                        autopas::Configuration(autopas::ContainerOption::verletClusterLists, 1.,
                                                               autopas::TraversalOption::c01CombinedSoA,
                                                               autopas::DataLayoutOption::soa,
                                                               autopas::Newton3Option::enabled)});

  fullSearchMPI.addEvidence(1);
  fullSearchMPI.tune();
  fullSearchMPI.addEvidence(0);
  fullSearchMPI.tune();
  fullSearchMPI.addEvidence(2);
  fullSearchMPI.tune();

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1.2,
                                   autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::aos,
                                   autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());

  fullSearchMPI.reset();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::directSum, 1.,
                                   autopas::TraversalOption::directSumTraversal,
                                   autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());
}

TEST_F(FullSearchMPITest, testInvalidConfigs) {
  autopas::FullSearchMPI fullSearchMPI({autopas::Configuration(autopas::ContainerOption::directSum, 1.,
                                                               autopas::TraversalOption::directSumTraversal,
                                                               autopas::DataLayoutOption::soa,
                                                               autopas::Newton3Option::enabled),
                                        autopas::Configuration(autopas::ContainerOption::linkedCells, 1.2,
                                                               autopas::TraversalOption::c08,
                                                               autopas::DataLayoutOption::aos,
                                                               autopas::Newton3Option::enabled),
                                        autopas::Configuration(autopas::ContainerOption::verletClusterLists, 1.,
                                                               autopas::TraversalOption::c01CombinedSoA,
                                                               autopas::DataLayoutOption::soa,
                                                               autopas::Newton3Option::enabled)});

  fullSearchMPI.addEvidence(1);
  fullSearchMPI.tune();
  fullSearchMPI.addEvidence(0);
  fullSearchMPI.tune();
  fullSearchMPI.addEvidence(2);
  fullSearchMPI.tune();

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1.2,
                                   autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::aos,
                                   autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());

  fullSearchMPI.tune(true);

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::directSum, 1.,
                                   autopas::TraversalOption::directSumTraversal,
                                   autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());

  fullSearchMPI.tune(true);

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::verletClusterLists, 1.,
                                   autopas::TraversalOption::c01CombinedSoA,
                                   autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::enabled),
            fullSearchMPI.getCurrentConfiguration());
}
