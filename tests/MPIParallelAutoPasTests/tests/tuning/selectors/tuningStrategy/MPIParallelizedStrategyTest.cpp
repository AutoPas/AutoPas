/**
 * @file MPIParallelizedStrategyTest.cpp
 * @author W. Thieme
 * @date 11.06.2020
 */

#include "MPIParallelizedStrategyTest.h"

#include <gmock/gmock-matchers.h>

#include "autopas/tuning/tuningStrategy/MPIParallelizedStrategy.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"

void MPIParallelizedStrategyTest::testBucketDistribution(
    const std::array<double, numRanksExpected> &particleDependentBinDensityStdDevs,
    const std::array<size_t, numRanksExpected> &expectedNumLocalConfigs,
    const std::set<autopas::Configuration> &searchSpace) {
  // Get MPI info
  int rank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
  int numRanks{};
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);

  // Setup Sanity checks
  ASSERT_EQ(numRanks, numRanksExpected) << "This test expects there to be three communicating MPI ranks!";
  ASSERT_LE(searchSpace.size(), std::accumulate(expectedNumLocalConfigs.begin(), expectedNumLocalConfigs.end(), 0ul))
      << "The sum of expected local configurations has to match the search space size times the number of buckets.";
  ASSERT_GT(searchSpace.size(), numRanks) << "There should be more configurations than ranks, otherwise the test can't "
                                             "distribute them in the desired way.";

  constexpr double mpiTuningMaxDifferenceForBucket{0.3};
  constexpr double mpiTuningWeightForMaxDensity{0.0};
  const autopas::AutoPas_MPI_Comm autopasMpiCommunicator = AUTOPAS_MPI_COMM_WORLD;

  auto mpiDivideAndConquerStrategy = autopas::MPIParallelizedStrategy(
      autopas::MPIParallelizedStrategy::createFallBackConfiguration(searchSpace,
                                                                    autopas::InteractionTypeOption::pairwise),
      autopasMpiCommunicator, mpiTuningMaxDifferenceForBucket, mpiTuningWeightForMaxDensity);

  // Begin with the full search space as the configQueue
  std::vector<autopas::Configuration> myConfigQueue{searchSpace.begin(), searchSpace.end()};

  // Send domain similarity statistics
  mpiDivideAndConquerStrategy.receiveDomainSimilarityStatistics(particleDependentBinDensityStdDevs[rank],
                                                                particleDependentBinDensityStdDevs[rank]);

  // Reset, which divides the configQueue between ranks which share a bucket
  mpiDivideAndConquerStrategy.reset(0, 0, myConfigQueue, {});

  // Check that all queues have the right size and the queues of each bucket add up to the search space
  EXPECT_EQ(myConfigQueue.size(), expectedNumLocalConfigs[rank])
      << "The local rank was not assigned the expected number of configurations.";

  // For every bucket check that the union of all configQueues adds up to the full search space
  const auto bucketCommunicator = mpiDivideAndConquerStrategy.getBucket();
  const auto globalConfs =
      autopas::utils::AutoPasConfigurationCommunicator::gatherConfigurations(bucketCommunicator, myConfigQueue, 0);
  int rankInBucket{};
  autopas::AutoPas_MPI_Comm_rank(bucketCommunicator, &rankInBucket);
  if (rankInBucket == 0) {
    EXPECT_GE(globalConfs.size(), searchSpace.size())
        << "Rank: " << rank << " RankInBucket: " << rankInBucket
        << " All rank's search spaces combined should be at least as many "
           "as the search space. It could be a bit more due to overlap.";
    const std::set<autopas::Configuration> uniqueConfigs(globalConfs.begin(), globalConfs.end());
    EXPECT_EQ(uniqueConfigs, searchSpace);
  }
}

/**
 * MPI divide and conquer for with one bucket.
 */
TEST_F(MPIParallelizedStrategyTest, testDistributeOneBucket) {
  const std::array<double, numRanksExpected> particleDependentBinDensityStdDevs{1., 1., 1.};
  const std::array<size_t, numRanksExpected> expectedNumLocalConfig{2, 2, 1};
  const std::set<autopas::Configuration> searchSpace{
      lc_c01_aos, lc_c04_aos, lc_c08_aos, lc_c01_soa, lc_c04_soa,
  };
  testBucketDistribution(particleDependentBinDensityStdDevs, expectedNumLocalConfig, searchSpace);
}

/**
 * MPI divide and conquer for with two buckets.
 */
TEST_F(MPIParallelizedStrategyTest, testDistributeTwoBucket) {
  const std::array<double, numRanksExpected> particleDependentBinDensityStdDevs{1., 10., 10.};
  const std::array<size_t, numRanksExpected> expectedNumLocalConfig{5, 3, 2};
  const std::set<autopas::Configuration> searchSpace{
      lc_c01_aos, lc_c04_aos, lc_c08_aos, lc_c01_soa, lc_c04_soa,
  };
  testBucketDistribution(particleDependentBinDensityStdDevs, expectedNumLocalConfig, searchSpace);
}
