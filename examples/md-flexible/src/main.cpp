/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "autopas/utils/WrapMPI.h"
#include "Simulation.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  autopas::AutoPas_MPI_Init(&argc, &argv);
  Simulation simulation(3, argc, argv);
  simulation.run();

  autopas::AutoPas_MPI_Finalize();
  return EXIT_SUCCESS;
}
