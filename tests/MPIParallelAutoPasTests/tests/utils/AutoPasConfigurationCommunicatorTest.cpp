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

// Test if the search space does get reduced.
TEST_F(AutoPasConfigurationCommunicatorTest, testDistributeConfigurationsFiniteCellSizeFactors) {
  std::set<ContainerOption> containerOptions{ContainerOption::verletClusterLists, ContainerOption::linkedCells};
  NumberSetFinite<double> cellSizeFactors{0.9, 1.0, 1.1};
  std::set<TraversalOption> traversalOptions{TraversalOption::vcl_cluster_iteration, TraversalOption::lc_sliced};
  std::set<LoadEstimatorOption> loadEstimatorOptions{LoadEstimatorOption::none,
                                                     LoadEstimatorOption::squaredParticlesPerCell};
  std::set<DataLayoutOption> dataLayoutOptions{DataLayoutOption::aos, DataLayoutOption::soa};
  std::set<Newton3Option> newton3Options{Newton3Option::enabled, Newton3Option::disabled};
  std::set<VectorizationPatternOption> vecPatternOptions{VectorizationPatternOption::p1xVec,
                                                         VectorizationPatternOption::pVecDiv2x2};
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  int totalNumConfigsBefore =
      getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions, dataLayoutOptions,
                         newton3Options, InteractionTypeOption::pairwise, vecPatternOptions);
  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions, dataLayoutOptions,
                           newton3Options, InteractionTypeOption::pairwise, vecPatternOptions, rank, commSize);
  int totalNumConfigsAfter =
      getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions, dataLayoutOptions,
                         newton3Options, InteractionTypeOption::pairwise, vecPatternOptions);

  // If true, each rank should have several configurations left.
  if (commSize <= totalNumConfigsBefore) {
    // No upper bound can properly be tested, because configurations are converted to sets of options, inflating the
    // search space.
    EXPECT_GE(totalNumConfigsAfter, totalNumConfigsBefore / commSize);
  } else {
    EXPECT_EQ(totalNumConfigsAfter, 1);
  }
}

// tests the precise outcomes of each rank for a fictional commSize of 4
// The distribution has been computed manually.
TEST_F(AutoPasConfigurationCommunicatorTest, testDistributeConfigurationsFiniteCellSizeFactorsFakeMPI) {
  std::set<ContainerOption> containerOptions{ContainerOption::verletClusterLists, ContainerOption::linkedCells};
  NumberSetFinite<double> cellSizeFactors{0.9, 1.0, 1.1};
  std::set<TraversalOption> traversalOptions{TraversalOption::vcl_cluster_iteration, TraversalOption::lc_sliced};
  std::set<LoadEstimatorOption> loadEstimatorOptions{LoadEstimatorOption::none};
  std::set<DataLayoutOption> dataLayoutOptions{DataLayoutOption::aos, DataLayoutOption::soa};
  std::set<Newton3Option> newton3Options{Newton3Option::enabled, Newton3Option::disabled};
  std::set<VectorizationPatternOption> vecPatternOptions{VectorizationPatternOption::p2xVecDiv2,
                                                         VectorizationPatternOption::pVecx1};

  // Rank 0
  auto containersTmp = std::set<ContainerOption>(containerOptions);
  auto cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  auto traversalsTmp = std::set<TraversalOption>(traversalOptions);
  auto loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  auto dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  auto newton3Tmp = std::set<Newton3Option>(newton3Options);
  auto vecPatternsTmp = std::set<VectorizationPatternOption>(vecPatternOptions);
  auto firstAndSecondCellSizes = std::set<double>{0.9, 1.0};
  auto secondAndThirdCellSizes = std::set<double>{1.0, 1.1};

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorTmp, dataLayoutTmp,
                           newton3Tmp, InteractionTypeOption::pairwise, vecPatternsTmp, 0, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::linkedCells});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), firstAndSecondCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::lc_sliced});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);
  EXPECT_EQ(vecPatternsTmp, vecPatternOptions);

  // Rank 1
  containersTmp = std::set<ContainerOption>(containerOptions);
  cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  traversalsTmp = std::set<TraversalOption>(traversalOptions);
  loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  newton3Tmp = std::set<Newton3Option>(newton3Options);
  vecPatternsTmp = std::set<VectorizationPatternOption>(vecPatternOptions);

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorTmp, dataLayoutTmp,
                           newton3Tmp, InteractionTypeOption::pairwise, vecPatternsTmp, 1, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::linkedCells});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), secondAndThirdCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::lc_sliced});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);
  EXPECT_EQ(vecPatternsTmp, vecPatternOptions);

  // Rank 2
  containersTmp = std::set<ContainerOption>(containerOptions);
  cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  traversalsTmp = std::set<TraversalOption>(traversalOptions);
  loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  newton3Tmp = std::set<Newton3Option>(newton3Options);
  vecPatternsTmp = std::set<VectorizationPatternOption>(vecPatternOptions);

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorOptions, dataLayoutTmp,
                           newton3Tmp, InteractionTypeOption::pairwise, vecPatternsTmp, 2, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::verletClusterLists});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), firstAndSecondCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::vcl_cluster_iteration});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);
  EXPECT_EQ(vecPatternsTmp, vecPatternOptions);

  // Rank 3
  containersTmp = std::set<ContainerOption>(containerOptions);
  cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  traversalsTmp = std::set<TraversalOption>(traversalOptions);
  loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  newton3Tmp = std::set<Newton3Option>(newton3Options);
  vecPatternsTmp = std::set<VectorizationPatternOption>(vecPatternOptions);

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorTmp, dataLayoutTmp,
                           newton3Tmp, InteractionTypeOption::pairwise, vecPatternsTmp, 3, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::verletClusterLists});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), secondAndThirdCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::vcl_cluster_iteration});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);
  EXPECT_EQ(vecPatternsTmp, vecPatternOptions);
}

// Test if CSFs are distributed if only one configuration exists for several ranks.
TEST_F(AutoPasConfigurationCommunicatorTest, testDistributeConfigurationsInfiniteCellSizeFactors) {
  std::set<ContainerOption> containerOptions{ContainerOption::verletClusterLists};
  NumberInterval<double> cellSizeFactors{0.8, 1.2};
  std::set<TraversalOption> traversalOptions{TraversalOption::vcl_cluster_iteration};
  std::set<LoadEstimatorOption> loadEstimatorOptions{LoadEstimatorOption::squaredParticlesPerCell};
  std::set<DataLayoutOption> dataLayoutOptions{DataLayoutOption::aos};
  std::set<Newton3Option> newton3Options{Newton3Option::enabled};
  std::set<VectorizationPatternOption> vecPatternOptions{VectorizationPatternOption::p1xVec};
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions, dataLayoutOptions,
                           newton3Options, InteractionTypeOption::pairwise, vecPatternOptions, rank, commSize);

  // Distribution should never return an empty search space.
  EXPECT_FALSE(containerOptions.empty() or cellSizeFactors.isEmpty() or traversalOptions.empty() or
               dataLayoutOptions.empty() or newton3Options.empty() or vecPatternOptions.empty());
  double error = 0.001;
  // Test even distribution.
  // Example of rank = 0 for 4 ranks in total:
  // Bounds: 0.8 + (0.4 / 4) * 0 = 0.8 and 0.8 + (0.4 / 4) * 1 = 0.9
  EXPECT_GE(cellSizeFactors.getMin(), 0.8 + (0.4 / commSize) * rank - error);
  EXPECT_LE(cellSizeFactors.getMin(), 0.8 + (0.4 / commSize) * rank + error);
  EXPECT_GE(cellSizeFactors.getMax(), 0.8 + (0.4 / commSize) * (rank + 1) - error);
  EXPECT_LE(cellSizeFactors.getMax(), 0.8 + (0.4 / commSize) * (rank + 1) + error);
}

// Test where the number of configurations equals the number of ranks.
TEST_F(AutoPasConfigurationCommunicatorTest, testDistributeOneConfigPerRank) {
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  std::set<ContainerOption> oneContainer{ContainerOption::linkedCells};
  std::set<double> cellSizeSet;
  for (int i = 0; i < commSize; ++i) {
    cellSizeSet.emplace(1.0 + i / 100.0);
  }
  autopas::NumberSetFinite rankManyCellSizes(cellSizeSet);
  std::set<TraversalOption> oneTraversal{TraversalOption::lc_c08};
  std::set<LoadEstimatorOption> oneLoadEstimator{LoadEstimatorOption::none};
  std::set<DataLayoutOption> oneDataLayout{DataLayoutOption::aos};
  std::set<Newton3Option> oneNewton3{Newton3Option::disabled};
  std::set<VectorizationPatternOption> oneVecPattern{VectorizationPatternOption::pVecx1};

  distributeConfigurations(oneContainer, rankManyCellSizes, oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3,
                           InteractionTypeOption::pairwise, oneVecPattern, rank, commSize);
  size_t size = getSearchSpaceSize(oneContainer, rankManyCellSizes, oneTraversal, oneLoadEstimator, oneDataLayout,
                                   oneNewton3, InteractionTypeOption::pairwise, oneVecPattern);

  EXPECT_EQ(size, 1);
  double error = 0.001;
  // Assumes that the CSFs are ordered increasingly in cellSizeSet.
  EXPECT_GE(rankManyCellSizes.getAll(), std::set<double>{1.0 + rank / 100.0 - error});
  EXPECT_LE(rankManyCellSizes.getAll(), std::set<double>{1.0 + rank / 100.0 + error});
}

// Example search space.
TEST_F(AutoPasConfigurationCommunicatorTest, testGetSearchSpaceSizeValid) {
  std::set<ContainerOption> threeContainers{ContainerOption::linkedCells, ContainerOption::verletClusterLists,
                                            ContainerOption::directSum};
  autopas::NumberSetFinite<double> twoCellSizes{1, 1.2};
  std::set<TraversalOption> threeTraversals{TraversalOption::lc_c08, TraversalOption::vcl_c06,
                                            TraversalOption::ds_sequential};
  std::set<LoadEstimatorOption> twoLoadEstimators{LoadEstimatorOption::none,
                                                  LoadEstimatorOption::squaredParticlesPerCell};
  std::set<DataLayoutOption> oneDataLayout{DataLayoutOption::aos};
  std::set<Newton3Option> oneNewton3{Newton3Option::disabled};
  std::set<VectorizationPatternOption> threeVecPatterns{VectorizationPatternOption::p2xVecDiv2,
                                                        VectorizationPatternOption::pVecDiv2x2,
                                                        VectorizationPatternOption::pVecx1};

  size_t size = getSearchSpaceSize(threeContainers, twoCellSizes, threeTraversals, twoLoadEstimators, oneDataLayout,
                                   oneNewton3, InteractionTypeOption::pairwise, threeVecPatterns);

  // There are 108 configurations in the Cartesian product, but only 18 of them are valid.
  EXPECT_EQ(size, 18);
}

// Example search space without valid configurations.
TEST_F(AutoPasConfigurationCommunicatorTest, testGetSearchSpaceSizeInvalid) {
  std::set<ContainerOption> twoContainers{ContainerOption::linkedCells, ContainerOption::verletListsCells};
  autopas::NumberSetFinite<double> twoCellSizes{1, 1.2};
  std::set<TraversalOption> twoTraversals{TraversalOption::ds_sequential, TraversalOption::vcl_cluster_iteration};
  std::set<LoadEstimatorOption> oneLoadEstimators{LoadEstimatorOption::neighborListLength};
  std::set<DataLayoutOption> oneDataLayout{DataLayoutOption::aos};
  std::set<Newton3Option> oneNewton3{Newton3Option::disabled};
  std::set<VectorizationPatternOption> oneVecPattern{VectorizationPatternOption::p1xVec};

  size_t size = getSearchSpaceSize(twoContainers, twoCellSizes, twoTraversals, oneLoadEstimators, oneDataLayout,
                                   oneNewton3, InteractionTypeOption::pairwise, oneVecPattern);

  // There are 8 configurations in the Cartesian product, but none are valid.
  EXPECT_EQ(size, 0);
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