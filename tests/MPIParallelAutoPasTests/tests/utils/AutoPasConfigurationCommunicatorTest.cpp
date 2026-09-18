/**
 * @file AutoPasConfigurationCommunicatorTest.h
 * @author W. Thieme
 * @date 05.06.2020
 */

#include "AutoPasConfigurationCommunicatorTest.h"

using namespace autopas::utils::AutoPasConfigurationCommunicator;
using namespace autopas;

// Test if serializing and deserializing again works as expected.
TEST_F(AutoPasConfigurationCommunicatorTest, testSerializeAndDeserialize) {
  Configuration config = Configuration(ContainerOption::directSum, 1.2, TraversalOption::lc_sliced,
                                       LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::disabled,
                                       InteractionTypeOption::pairwise, VectorizationPatternOption::p1xVec);
  Configuration passedConfig = deserializeConfiguration(serializeConfiguration(config));
  EXPECT_EQ(passedConfig, config);
}

// Test if serializing and deserializing a vector of configurations works as expected.
TEST_F(AutoPasConfigurationCommunicatorTest, testSerializeAndDeserializeVector) {
  const std::vector<autopas::Configuration> configurations = {
      autopas::Configuration{autopas::ContainerOption::octree, 1., autopas::TraversalOption::ot_c18,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                             autopas::Newton3Option::disabled, InteractionTypeOption::pairwise,
                             VectorizationPatternOption::p1xVec},
      autopas::Configuration{autopas::ContainerOption::verletClusterLists, 1., autopas::TraversalOption::vcl_c06,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, InteractionTypeOption::pairwise,
                             VectorizationPatternOption::pVecDiv2x2},
      autopas::Configuration{autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced_balanced,
                             autopas::LoadEstimatorOption::squaredParticlesPerCell, autopas::DataLayoutOption::aos,
                             autopas::Newton3Option::enabled, InteractionTypeOption::pairwise,
                             VectorizationPatternOption::pVecx1},
  };
  const auto serializedConfigs = serializeConfigurations(configurations);
  const auto passedConfig = deserializeConfigurations(serializedConfigs);
  EXPECT_EQ(passedConfig, configurations);
}

// Test if the optimization distributes the configuration with the lowest provided time.
TEST_F(AutoPasConfigurationCommunicatorTest, testOptimizeConfiguration) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Configuration config =
      Configuration(ContainerOption::directSum, 1 + rank, TraversalOption::lc_sliced,
                    LoadEstimatorOption::neighborListLength, DataLayoutOption::aos, Newton3Option::enabled,
                    InteractionTypeOption::pairwise, VectorizationPatternOption::p1xVec);
  // provide rank as the time for the config.
  Configuration optimized = findGloballyBestConfiguration(MPI_COMM_WORLD, config, rank);

  // CSF should be 1, because rank 0 provided the lowest time.
  EXPECT_EQ(optimized,
            Configuration(ContainerOption::directSum, 1, TraversalOption::lc_sliced,
                          LoadEstimatorOption::neighborListLength, DataLayoutOption::aos, Newton3Option::enabled,
                          InteractionTypeOption::pairwise, VectorizationPatternOption::p1xVec));
}

TEST_F(AutoPasConfigurationCommunicatorTest, testGatherConfigs) {
  int rank{};
  AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
  constexpr int numRanksExpected = 3;
  int numRanks{};
  AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);

  ASSERT_EQ(numRanks, numRanksExpected) << "This test expects there to be three communicating MPI ranks!";

  const std::vector<Configuration> expectedConfigurations{
      autopas::Configuration{autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                             autopas::Newton3Option::disabled, InteractionTypeOption::pairwise,
                             VectorizationPatternOption::p1xVec},
      autopas::Configuration{autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c04,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                             autopas::Newton3Option::disabled, InteractionTypeOption::pairwise,
                             VectorizationPatternOption::p1xVec},
      autopas::Configuration{autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                             autopas::Newton3Option::disabled, InteractionTypeOption::pairwise,
                             VectorizationPatternOption::p1xVec},
  };

  const auto localConf = [&]() -> std::vector<Configuration> {
    switch (rank) {
      case 0: {
        return {
            expectedConfigurations[0],
        };
      }
      case 1: {
        return {
            expectedConfigurations[1],
        };
      }
      case 2: {
        return {
            expectedConfigurations[2],
        };
      }
      default: {
        // should never happen.
        return {};
      };
    }
  }();

  const auto gatheredConfigurations =
      autopas::utils::AutoPasConfigurationCommunicator::gatherConfigurations(AUTOPAS_MPI_COMM_WORLD, localConf, 0);

  if (rank == 0) {
    EXPECT_EQ(expectedConfigurations, gatheredConfigurations);
  }
}