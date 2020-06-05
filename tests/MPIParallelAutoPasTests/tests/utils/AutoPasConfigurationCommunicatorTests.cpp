/**
 * @file AutoPasConfigurationCommunicator.h
 * @author W. Thieme
 * @date 05.06.2020
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include "autopas/utils/AutoPasConfigurationCommunicator.h"

using namespace autopas::AutoPasConfigurationCommunicator;
using namespace autopas;

TEST(AutoPasConfigurationCommunicatorTest, testSerializeAndDeserialize) {
  Configuration config = Configuration(ContainerOption::directSum, 1.2, TraversalOption::sliced, DataLayoutOption::cuda,
                                       Newton3Option::disabled);
  Configuration passedConfig = deserializeConfiguration(serializeConfiguration(config));
  EXPECT_EQ(passedConfig, config);
}

TEST(AutoPasConfigurationCommunicatorTest, testOptimizeConfiguration) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Configuration config = Configuration(ContainerOption::directSum, 1 + rank, TraversalOption::sliced,
                                       DataLayoutOption::aos, Newton3Option::enabled);
  Configuration optimized = optimizeConfiguration(MPI_COMM_WORLD, config, rank);

  EXPECT_EQ(optimized, Configuration(ContainerOption::directSum, 1, TraversalOption::sliced, DataLayoutOption::aos,
                                     Newton3Option::enabled));
}

TEST(AutoPasConfigurationCommunicatorTest, testDistributeConfigurationsFiniteCellSizeFactors) {
  std::set<ContainerOption> containerOptions{ContainerOption::verletClusterLists, ContainerOption::linkedCells};
  NumberSetFinite<double> cellSizeFactors{0.9, 1.0, 1.1};
  std::set<TraversalOption> traversalOptions{TraversalOption::verletClusters, TraversalOption::sliced};
  std::set<DataLayoutOption> dataLayoutOptions{DataLayoutOption::aos, DataLayoutOption::soa};
  std::set<Newton3Option> newton3Options{Newton3Option::enabled, Newton3Option::disabled};
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions, newton3Options,
                           MPI_COMM_WORLD);

  EXPECT_FALSE(containerOptions.empty() or cellSizeFactors.isEmpty() or traversalOptions.empty() or
               dataLayoutOptions.empty() or newton3Options.empty());

  // multiplying all options together gives 48 possible configurations
  EXPECT_GE(containerOptions.size() * cellSizeFactors.size() * traversalOptions.size() * dataLayoutOptions.size()
            * newton3Options.size(), 48/commSize);

  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions, newton3Options,
                           MPI_COMM_SELF);

  EXPECT_FALSE(containerOptions.empty() or cellSizeFactors.isEmpty() or traversalOptions.empty() or
               dataLayoutOptions.empty() or newton3Options.empty());
}

TEST(AutoPasConfigurationCommunicatorTest, testDistributeConfigurationsInfiniteCellSizeFactors) {
  std::set<ContainerOption> containerOptions{ContainerOption::verletClusterLists};
  NumberInterval<double> cellSizeFactors{0.8, 1.2};
  std::set<TraversalOption> traversalOptions{TraversalOption::verletClusters};
  std::set<DataLayoutOption> dataLayoutOptions{DataLayoutOption::aos};
  std::set<Newton3Option> newton3Options{Newton3Option::enabled};
  int rank, commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  distributeConfigurations(containerOptions, cellSizeFactors, traversalOptions, dataLayoutOptions, newton3Options,
                           MPI_COMM_WORLD);

  EXPECT_FALSE(containerOptions.empty() or cellSizeFactors.isEmpty() or traversalOptions.empty() or
               dataLayoutOptions.empty() or newton3Options.empty());
  double error = 0.001;
  EXPECT_GE(cellSizeFactors.getMin(), 0.8 + (0.4 / commSize) * rank - error);
  EXPECT_LE(cellSizeFactors.getMin(), 0.8 + (0.4 / commSize) * rank + error);
  EXPECT_GE(cellSizeFactors.getMax(), 0.8 + (0.4 / commSize) * (rank + 1) - error);
  EXPECT_LE(cellSizeFactors.getMax(), 0.8 + (0.4 / commSize) * (rank + 1) + error);
}
