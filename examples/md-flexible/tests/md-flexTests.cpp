/**
 * @file md-flexTests.cpp
 * @author F. Gratl
 * @date 06.11.19
 */

#include <gtest/gtest.h>

#include "autopas/utils/WrapMPI.h"

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  autopas::AutoPas_MPI_Init(&argc, &argv);
  return RUN_ALL_TESTS();
  autopas::AutoPas_MPI_Finalize();
}
