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
                                       LoadEstimatorOption::none, DataLayoutOption::cuda, Newton3Option::disabled);
  Configuration passedConfig = deserializeConfiguration(serializeConfiguration(config));
  EXPECT_EQ(passedConfig, config);
}

// Test if the optimization distributes the configuration with the lowest provided time.
TEST_F(AutoPasConfigurationCommunicatorTest, testOptimizeConfiguration) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Configuration config =
      Configuration(ContainerOption::directSum, 1 + rank, TraversalOption::lc_sliced,
                    LoadEstimatorOption::neighborListLength, DataLayoutOption::aos, Newton3Option::enabled);
  // provide rank as the time for the config.
  Configuration optimized = optimizeConfiguration(MPI_COMM_WORLD, config, rank);

  // CSF should be 1, because rank 0 provided the lowest time.
  EXPECT_EQ(optimized,
            Configuration(ContainerOption::directSum, 1, TraversalOption::lc_sliced,
                          LoadEstimatorOption::neighborListLength, DataLayoutOption::aos, Newton3Option::enabled));
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
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  int totalNumConfigsBefore = getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions,
                                                 loadEstimatorOptions, dataLayoutOptions, newton3Options);
  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions, dataLayoutOptions,
                           newton3Options, rank, commSize);
  int totalNumConfigsAfter = getSearchSpaceSize(containerOptions, cellSizeFactors, traversalOptions,
                                                loadEstimatorOptions, dataLayoutOptions, newton3Options);

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

  // Rank 0
  auto containersTmp = std::set<ContainerOption>(containerOptions);
  auto cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  auto traversalsTmp = std::set<TraversalOption>(traversalOptions);
  auto loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  auto dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  auto newton3Tmp = std::set<Newton3Option>(newton3Options);
  auto firstAndSecondCellSizes = std::set<double>{0.9, 1.0};
  auto secondAndThirdCellSizes = std::set<double>{1.0, 1.1};

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorTmp, dataLayoutTmp,
                           newton3Tmp, 0, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::linkedCells});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), firstAndSecondCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::lc_sliced});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);

  // Rank 1
  containersTmp = std::set<ContainerOption>(containerOptions);
  cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  traversalsTmp = std::set<TraversalOption>(traversalOptions);
  loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  newton3Tmp = std::set<Newton3Option>(newton3Options);

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorTmp, dataLayoutTmp,
                           newton3Tmp, 1, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::linkedCells});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), secondAndThirdCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::lc_sliced});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);

  // Rank 2
  containersTmp = std::set<ContainerOption>(containerOptions);
  cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  traversalsTmp = std::set<TraversalOption>(traversalOptions);
  loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  newton3Tmp = std::set<Newton3Option>(newton3Options);

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorOptions, dataLayoutTmp,
                           newton3Tmp, 2, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::verletClusterLists});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), firstAndSecondCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::vcl_cluster_iteration});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);

  // Rank 3
  containersTmp = std::set<ContainerOption>(containerOptions);
  cellSizeFactorsTmp = NumberSetFinite<double>(cellSizeFactors);
  traversalsTmp = std::set<TraversalOption>(traversalOptions);
  loadEstimatorTmp = std::set<LoadEstimatorOption>(loadEstimatorOptions);
  dataLayoutTmp = std::set<DataLayoutOption>(dataLayoutOptions);
  newton3Tmp = std::set<Newton3Option>(newton3Options);

  distributeConfigurations(containersTmp, cellSizeFactorsTmp, traversalsTmp, loadEstimatorTmp, dataLayoutTmp,
                           newton3Tmp, 3, 4);
  EXPECT_EQ(containersTmp, std::set<ContainerOption>{ContainerOption::verletClusterLists});
  EXPECT_EQ(cellSizeFactorsTmp.getAll(), secondAndThirdCellSizes);
  EXPECT_EQ(traversalsTmp, std::set<TraversalOption>{TraversalOption::vcl_cluster_iteration});
  EXPECT_EQ(loadEstimatorTmp, std::set<LoadEstimatorOption>{LoadEstimatorOption::none});
  EXPECT_EQ(dataLayoutTmp, dataLayoutOptions);
  EXPECT_EQ(newton3Tmp, newton3Options);
}

// Test if CSFs are distributed if only one configuration exists for several ranks.
TEST_F(AutoPasConfigurationCommunicatorTest, testDistributeConfigurationsInfiniteCellSizeFactors) {
  std::set<ContainerOption> containerOptions{ContainerOption::verletClusterLists};
  NumberInterval<double> cellSizeFactors{0.8, 1.2};
  std::set<TraversalOption> traversalOptions{TraversalOption::vcl_cluster_iteration};
  std::set<LoadEstimatorOption> loadEstimatorOptions{LoadEstimatorOption::squaredParticlesPerCell};
  std::set<DataLayoutOption> dataLayoutOptions{DataLayoutOption::aos};
  std::set<Newton3Option> newton3Options{Newton3Option::enabled};
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, loadEstimatorOptions, dataLayoutOptions,
                           newton3Options, rank, commSize);

  // Distribution should never return an empty search space.
  EXPECT_FALSE(containerOptions.empty() or cellSizeFactors.isEmpty() or traversalOptions.empty() or
               dataLayoutOptions.empty() or newton3Options.empty());
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

  distributeConfigurations(oneContainer, rankManyCellSizes, oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3,
                           rank, commSize);
  size_t size =
      getSearchSpaceSize(oneContainer, rankManyCellSizes, oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3);

  EXPECT_EQ(size, 1);
  double error = 0.001;
  // Assumes that the CSFs are ordered increasingly in cellSizeSet.
  EXPECT_GE(rankManyCellSizes.getAll(), std::set<double>{1.0 + rank / 100.0 - error});
  EXPECT_LE(rankManyCellSizes.getAll(), std::set<double>{1.0 + rank / 100.0 + error});
}

// Example search space.
TEST_F(AutoPasConfigurationCommunicatorTest, testGetSearchSpaceSizeValid) {
  std::set<ContainerOption> threeContainers{ContainerOption::linkedCells, ContainerOption::verletClusterCells,
                                            ContainerOption::directSum};
  autopas::NumberSetFinite<double> twoCellSizes{1, 1.2};
  std::set<TraversalOption> threeTraversals{TraversalOption::lc_c08, TraversalOption::vcc_cluster_iteration_cuda,
                                            TraversalOption::ds_sequential};
  std::set<LoadEstimatorOption> twoLoadEstimators{LoadEstimatorOption::none,
                                                  LoadEstimatorOption::squaredParticlesPerCell};
  std::set<DataLayoutOption> oneDataLayout{DataLayoutOption::aos};
  std::set<Newton3Option> oneNewton3{Newton3Option::disabled};

  size_t size =
      getSearchSpaceSize(threeContainers, twoCellSizes, threeTraversals, twoLoadEstimators, oneDataLayout, oneNewton3);

  // There are 36 configurations in the Cartesian product, but only 6 of them are valid.
  EXPECT_EQ(size, 6);
}

// Example search space without valid configurations.
TEST_F(AutoPasConfigurationCommunicatorTest, testGetSearchSpaceSizeInvalid) {
  std::set<ContainerOption> twoContainers{ContainerOption::linkedCells, ContainerOption::verletListsCells};
  autopas::NumberSetFinite<double> twoCellSizes{1, 1.2};
  std::set<TraversalOption> twoTraversals{TraversalOption::ds_sequential, TraversalOption::vcl_cluster_iteration};
  std::set<LoadEstimatorOption> oneLoadEstimators{LoadEstimatorOption::neighborListLength};
  std::set<DataLayoutOption> oneDataLayout{DataLayoutOption::aos};
  std::set<Newton3Option> oneNewton3{Newton3Option::disabled};

  size_t size =
      getSearchSpaceSize(twoContainers, twoCellSizes, twoTraversals, oneLoadEstimators, oneDataLayout, oneNewton3);

  // There are 8 configurations in the Cartesian product, but none are valid.
  EXPECT_EQ(size, 0);
}
