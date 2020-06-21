/**
 * @file MPIParallelizedStrategyTest.cpp
 * @author W. Thieme
 * @date 11.06.2020
 */

#include "MPIParallelizedStrategyTest.h"

using namespace autopas;

/**
 * Simple setup for two configurations with different cellSizeFactors depending on rank
 * @param mpiParallelizedStrategy
 * @param rank
 */
void finiteCellSizeFactorsSetup(MPIParallelizedStrategy &mpiParallelizedStrategy) {
  do {
    long evidence = static_cast<long>(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor * 100.0);
    mpiParallelizedStrategy.addEvidence(evidence, 0);
  } while (mpiParallelizedStrategy.tune(false));
}

/**
 * Simple setup for configurations differing only in the min-max value of their cellSizeFactors
 * @param mpiParallelizedStrategy
 * @param rank
 * @return the smallest cellSizeFactor that was generated
 */
double infiniteCellSizeFactorSetup(MPIParallelizedStrategy &mpiParallelizedStrategy) {
  double smallestLocalCellSizeFactor = mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor;

  do {
    long evidence = static_cast<long>(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor * 100.0);

    mpiParallelizedStrategy.addEvidence(evidence, 0);

    if (smallestLocalCellSizeFactor > mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor) {
      smallestLocalCellSizeFactor = mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor;
    }
  } while (mpiParallelizedStrategy.tune(false));
  return smallestLocalCellSizeFactor;
}

/**
 * As it currently stands, these tests cannot work with AUTOPAS_MPI=OFF
 */

TEST_F(MPIParallelizedStrategyTest, testTuneFullSearch) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<FullSearch>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                               std::set<double>{1.0 + (2 * rank) / 10.0,
                                                                                   1.0 + (2 * rank + 1) / 10.0},
                                               std::set<TraversalOption>{TraversalOption::sliced},
                                               std::set<LoadEstimatorOption>{LoadEstimatorOption::none},
                                               std::set<DataLayoutOption>{DataLayoutOption::aos},
                                               std::set<Newton3Option>{Newton3Option::enabled});

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy);

  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::sliced, LoadEstimatorOption::none,DataLayoutOption::aos,
                          Newton3Option::enabled));
}

TEST_F(MPIParallelizedStrategyTest, testTuneRandomSearchFiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<RandomSearch>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                 NumberSetFinite<double>{1.0 + (2 * rank) / 10.0,
                                                                         1.0 + (2 * rank + 1) / 10.0},
                                                 std::set<TraversalOption>{TraversalOption::sliced},
                                                 std::set<LoadEstimatorOption>{LoadEstimatorOption::none},
                                                 std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                 std::set<Newton3Option>{Newton3Option::enabled},
                                                 2);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy);

  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, 1.1);
}

TEST_F(MPIParallelizedStrategyTest, testTuneRandomSearchInfiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<RandomSearch>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                 NumberInterval<double>{1.0 + (2 * rank) / 10.0,
                                                                        1.0 + (2 * rank + 1) / 10.0},
                                                 std::set<TraversalOption>{TraversalOption::sliced},
                                                 std::set<LoadEstimatorOption>{LoadEstimatorOption::none},
                                                 std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                 std::set<Newton3Option>{Newton3Option::enabled},
                                                 2);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  auto smallestLocalCellSizeFactor = infiniteCellSizeFactorSetup(mpiParallelizedStrategy);

  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, smallestLocalCellSizeFactor);
}

TEST_F(MPIParallelizedStrategyTest, testTuneActiveHarmonyFiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<ActiveHarmony>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                  NumberSetFinite<double>{1.0 + (2 * rank) / 10.0,
                                                                          1.0 + (2 * rank + 1) / 10.0},
                                                  std::set<TraversalOption>{TraversalOption::sliced},
                                                  std::set<LoadEstimatorOption>{LoadEstimatorOption::none},
                                                  std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                  std::set<Newton3Option>{Newton3Option::enabled});

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy);

  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::sliced, LoadEstimatorOption::none, DataLayoutOption::aos,
                          Newton3Option::enabled));
}

TEST_F(MPIParallelizedStrategyTest, testTuneBayesianClusterSearchInfiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<BayesianClusterSearch>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                          NumberSetFinite<double>{1.0 + (2 * rank) / 10.0,
                                                                                  1.0 + (2 * rank + 1) / 10.0},
                                                          std::set<TraversalOption>{TraversalOption::sliced},
                                                          std::set<LoadEstimatorOption>{LoadEstimatorOption::none},
                                                          std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                          std::set<Newton3Option>{Newton3Option::enabled});

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  auto smallestCellSizeFactor = infiniteCellSizeFactorSetup(mpiParallelizedStrategy);

  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, smallestCellSizeFactor);
}