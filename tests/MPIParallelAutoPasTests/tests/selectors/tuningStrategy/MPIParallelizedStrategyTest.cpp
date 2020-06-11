/**
 * @file MPIParallelizedStrategyTest.cpp
 * @author W. Thieme
 * @date 11.06.2020
 */

#include "MPIParallelizedStrategyTest.h"

using namespace autopas;

/**
 * Simple setup for two configurations witch different cellSizeFactors depending on rank
 * @param mpiParallelizedStrategy
 * @param rank
 */
void finiteCellSizeFactorsSetup(MPIParallelizedStrategy &mpiParallelizedStrategy, int rank) {
  // complicated, because I cannot assume any order in configurations
  // for every configuration with cellSizeFactor with cellSizeFactor 1.0 + (2 * rank + x) / 10,
  // the evidence should be 2 * rank + x (with x in {0, 1})
  int evidenceOffset =
      static_cast<int>((mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor - 1.0) * 10.0) - (2 * rank);
  mpiParallelizedStrategy.addEvidence(2 * rank + evidenceOffset, 0);
  mpiParallelizedStrategy.tune(false);
  evidenceOffset = 1 - evidenceOffset;
  mpiParallelizedStrategy.addEvidence(2 * rank + evidenceOffset, 0);
  mpiParallelizedStrategy.tune(false);

  // now the local search spaces should be finished, so tune one last time for global tuning
  mpiParallelizedStrategy.tune(false);
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
                                               std::set<DataLayoutOption>{DataLayoutOption::aos},
                                               std::set<Newton3Option>{Newton3Option::enabled});

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy, rank);

  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::sliced, DataLayoutOption::aos,
                          Newton3Option::enabled));
}

TEST_F(MPIParallelizedStrategyTest, testTuneRandomSearchFiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<RandomSearch>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                 NumberSetFinite<double>{1.0 + (2 * rank) / 10.0,
                                                                         1.0 + (2 * rank + 1) / 10.0},
                                                 std::set<TraversalOption>{TraversalOption::sliced},
                                                 std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                 std::set<Newton3Option>{Newton3Option::enabled});

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy, rank);

  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::sliced, DataLayoutOption::aos,
                          Newton3Option::enabled));
}

TEST_F(MPIParallelizedStrategyTest, testTuneRandomSearchInfiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<RandomSearch>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                 NumberInterval<double>{1.0 + (2 * rank) / 10.0,
                                                                        1.0 + (2 * rank + 1) / 10.0},
                                                 std::set<TraversalOption>{TraversalOption::sliced},
                                                 std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                 std::set<Newton3Option>{Newton3Option::enabled},
                                                 2);

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);

  double smallestLocalCellSizeFactor = mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor;
  int evidence = static_cast<int>(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor * 100.0);

  mpiParallelizedStrategy.addEvidence(evidence, 0);
  mpiParallelizedStrategy.tune(false);

  if (smallestLocalCellSizeFactor > mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor) {
    smallestLocalCellSizeFactor = mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor;
  }
  evidence = static_cast<int>(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor * 100.0);
  mpiParallelizedStrategy.addEvidence(evidence, 0);
  mpiParallelizedStrategy.tune(false);

  // now the local search spaces should be finished, so tune one last time for global tuning
  mpiParallelizedStrategy.tune(false);

  EXPECT_LE(mpiParallelizedStrategy.getCurrentConfiguration().cellSizeFactor, smallestLocalCellSizeFactor);
}

TEST_F(MPIParallelizedStrategyTest, testTuneActiveHarmonyFiniteCellSizeFactors) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  auto strategy = std::make_unique<ActiveHarmony>(std::set<ContainerOption>{ContainerOption::linkedCells},
                                                  NumberSetFinite<double>{1.0 + (2 * rank) / 10.0,
                                                                          1.0 + (2 * rank + 1) / 10.0},
                                                  std::set<TraversalOption>{TraversalOption::sliced},
                                                  std::set<DataLayoutOption>{DataLayoutOption::aos},
                                                  std::set<Newton3Option>{Newton3Option::enabled});

  auto mpiParallelizedStrategy = MPIParallelizedStrategy(std::move(strategy), MPI_COMM_WORLD);
  finiteCellSizeFactorsSetup(mpiParallelizedStrategy, rank);

  EXPECT_EQ(mpiParallelizedStrategy.getCurrentConfiguration(),
            Configuration(ContainerOption::linkedCells, 1.0, TraversalOption::sliced, DataLayoutOption::aos,
                          Newton3Option::enabled));
}
