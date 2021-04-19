/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "simulations/MDFlexSingleNode.h"

/**
 * The main function for md-flexible.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv) {

	MDFlexSingleNode simulation(argc, argv);
	simulation.run();

  return EXIT_SUCCESS;
}
