/**
 * @file MPIParallelizedStrategyTest.cpp
 * @author W. Thieme
 * @date 11.06.2020
 */

#include "MPIParallelizedStrategyTest.h"

using namespace autopas;

/**
 * Simple setup for two configurations with different cellSizeFactors depending on rank.
 * The provided evidence increases strictly monotonically with the cellSizeFactor, so the configuration with the
 * smallest CSF wins.
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
 * Simple setup for configurations differing only in the min-max value of their cellSizeFactors.
 * The provided evidence increases strictly monotonically with the cellSizeFactor, so the configuration with the
 * smallest CSF wins.
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

TEST_F(MPIParallelizedStrategyTest, testTuneFullSearch) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto oneContainer = std::set<ContainerOption>{ContainerOption::linkedCells};
  const auto twoCellSizes = std::set<double>{1.0 + (2 * rank) / 10.0, 1.0 + (2 * rank + 1) / 10.0};
  const auto oneTraversal = std::set<TraversalOption>{TraversalOption::lc_sliced};
  const auto oneLoadEstimator = std::set<LoadEstimatorOption>{LoadEstimatorOption::none};
  const auto oneDataLayout = std::set<DataLayoutOption>{DataLayoutOption::aos};
  const auto oneNewton3 = std::set<Newton3Option>{Newton3Option::enabled};

  auto strategy = std::make_unique<FullSearch>(oneContainer, twoCellSizes, oneTraversal, oneLoadEstimator,
                                               oneDataLayout, oneNewton3);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD, oneContainer,
                                                         oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy);

  // CSF should be 1.0, because of how finiteCellSizeFactorsSetup() handles evidences.
  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::lc_sliced, LoadEstimatorOption::none,
                          DataLayoutOption::aos, Newton3Option::enabled));
}

TEST_F(MPIParallelizedStrategyTest, testTuneRandomSearchFiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto oneContainer = std::set<ContainerOption>{ContainerOption::linkedCells};
  const auto twoCellSizes = NumberSetFinite<double>{1.0 + (2 * rank) / 10.0, 1.0 + (2 * rank + 1) / 10.0};
  const auto oneTraversal = std::set<TraversalOption>{TraversalOption::lc_sliced};
  const auto oneLoadEstimator = std::set<LoadEstimatorOption>{LoadEstimatorOption::none};
  const auto oneDataLayout = std::set<DataLayoutOption>{DataLayoutOption::aos};
  const auto oneNewton3 = std::set<Newton3Option>{Newton3Option::enabled};

  auto strategy = std::make_unique<RandomSearch>(oneContainer, twoCellSizes, oneTraversal, oneLoadEstimator,
                                                 oneDataLayout, oneNewton3, 2);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD, oneContainer,
                                                         oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy);

  // Currently, RandomSearch ignores cellSizeFactors to determine how many configurations can be tested.
  // Thus, only one is tested despite maxEvidence=2, so we can only assume that the best configuration is at rank 0.
  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, 1.1);
}

TEST_F(MPIParallelizedStrategyTest, testTuneRandomSearchInfiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto oneContainer = std::set<ContainerOption>{ContainerOption::linkedCells};
  const auto twoCellSizes = NumberInterval<double>{1.0 + (2 * rank) / 10.0, 1.0 + (2 * rank + 1) / 10.0};
  const auto oneTraversal = std::set<TraversalOption>{TraversalOption::lc_sliced};
  const auto oneLoadEstimator = std::set<LoadEstimatorOption>{LoadEstimatorOption::none};
  const auto oneDataLayout = std::set<DataLayoutOption>{DataLayoutOption::aos};
  const auto oneNewton3 = std::set<Newton3Option>{Newton3Option::enabled};

  auto strategy = std::make_unique<RandomSearch>(oneContainer, twoCellSizes, oneTraversal, oneLoadEstimator,
                                                 oneDataLayout, oneNewton3, 2);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD, oneContainer,
                                                         oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3);
  auto smallestLocalCellSizeFactor = infiniteCellSizeFactorSetup(mpiParallelizedStrategy);

  // CSF should be the smallest, because of how finiteCellSizeFactorsSetup() handles evidences.
  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, smallestLocalCellSizeFactor);
}

TEST_F(MPIParallelizedStrategyTest, testTuneActiveHarmonyFiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto oneContainer = std::set<ContainerOption>{ContainerOption::linkedCells};
  const auto twoCellSizes = NumberSetFinite<double>{1.0 + (2 * rank) / 10.0, 1.0 + (2 * rank + 1) / 10.0};
  const auto oneTraversal = std::set<TraversalOption>{TraversalOption::lc_sliced};
  const auto oneLoadEstimator = std::set<LoadEstimatorOption>{LoadEstimatorOption::none};
  const auto oneDataLayout = std::set<DataLayoutOption>{DataLayoutOption::aos};
  const auto oneNewton3 = std::set<Newton3Option>{Newton3Option::enabled};

  auto strategy =
      std::make_unique<ActiveHarmony>(oneContainer, twoCellSizes, oneTraversal, oneLoadEstimator, oneDataLayout,
                                      oneNewton3, autopas::MPIStrategyOption::divideAndConquer, MPI_COMM_WORLD);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD, oneContainer,
                                                         oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy);

  // CSF should be 1.0, because of how finiteCellSizeFactorsSetup() handles evidences.
  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::lc_sliced, LoadEstimatorOption::none,
                          DataLayoutOption::aos, Newton3Option::enabled));
}

TEST_F(MPIParallelizedStrategyTest, testTuneBayesianSearchInfiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto oneContainer = std::set<ContainerOption>{ContainerOption::linkedCells};
  const auto twoCellSizes = NumberInterval<double>{1.0 + (2 * rank) / 10.0, 1.0 + (2 * rank + 1) / 10.0};
  const auto oneTraversal = std::set<TraversalOption>{TraversalOption::lc_sliced};
  const auto oneLoadEstimator = std::set<LoadEstimatorOption>{LoadEstimatorOption::none};
  const auto oneDataLayout = std::set<DataLayoutOption>{DataLayoutOption::aos};
  const auto oneNewton3 = std::set<Newton3Option>{Newton3Option::enabled};

  auto strategy = std::make_unique<BayesianSearch>(oneContainer, twoCellSizes, oneTraversal, oneLoadEstimator,
                                                   oneDataLayout, oneNewton3);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD, oneContainer,
                                                         oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3);
  auto smallestCellSizeFactor = infiniteCellSizeFactorSetup(mpiParallelizedStrategy);

  // BayesianSearch seems to not be absolutely optimal in this case.
  double error = 0.01;
  // CSF should be the smallest, because of how finiteCellSizeFactorsSetup() handles evidences.
  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, smallestCellSizeFactor + error);
}