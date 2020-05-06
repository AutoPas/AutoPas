/**
 * @file allTests.cpp
 * @author seckler
 * @date 18.01.18
 */

#include <gtest/gtest.h>
#include <mpi.h>

int main(int argc, char **argv) {
  int result = 0;

  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  MPI_Init(&argc, &argv);
  int worldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
  std::cout << "Testing rank " << worldRank << std::endl;

  // running only my tests
  result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
