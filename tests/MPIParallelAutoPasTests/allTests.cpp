/**
 * @file allTests.cpp
 * @author W. Thieme
 * @date 01.05.2020
 */

#include <gtest/gtest.h>

#include "autopas/utils/WrapMPI.h"

int main(int argc, char **argv) {
  int result = 0;

  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  autopas::AutoPas_MPI_Init(&argc, &argv);

  // CTest discovers the tests by running --gtest_list_tests under mpiexec. Only rank 0 may print the list, otherwise
  // every test is registered once per rank, and the interleaved output attaches tests to the wrong suites.
  int rank{0};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &rank);
  if (GTEST_FLAG_GET(list_tests) and rank != 0) {
    autopas::AutoPas_MPI_Finalize();
    return 0;
  }

  // running only my tests
  result = RUN_ALL_TESTS();

  autopas::AutoPas_MPI_Finalize();
  return result;
}
