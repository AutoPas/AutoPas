/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#if defined(AUTOPAS_INTERNODE_TUNING)
#include <mpi.h>
#endif

#include <iostream>

#include "simulations/MDFlexSingleNode.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {
	MDFlexSingleNode mdFlexSimulation;
	mdFlexSimulation.initialize(argc, argv);
	mdFlexSimulation.run();
	mdFlexSimulation.finalize(argc, argv);
  return EXIT_SUCCESS;
}
