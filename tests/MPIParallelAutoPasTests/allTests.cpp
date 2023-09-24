/**
 * @file allTests.cpp
 * @author W. Thieme
 * @date 01.05.2020
 */

#include <gtest/gtest.h>
#include <mpi.h>

#include "autopas/utils/WrapMPI.h"

int main(int argc, char **argv) {
  int result = 0;

  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  autopas::AutoPas_MPI_Init(&argc, &argv);

  // running only my tests
  result = RUN_ALL_TESTS();

  autopas::AutoPas_MPI_Finalize();
  return result;
}
