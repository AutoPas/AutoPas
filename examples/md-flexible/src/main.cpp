/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "simulations/MDFlexSingleNode.h"
#include "simulations/MDFlexMPI.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
  MDFlexMPI simulation(3, argc, argv);
  simulation.run();

  return EXIT_SUCCESS;
}
