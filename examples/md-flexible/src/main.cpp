/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "mpi.h"
#include "simulations/MDFlexMPI.h"
#include "simulations/MDFlexSingleRank.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MDFlexMPI simulation(3, argc, argv);
  simulation.run();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
