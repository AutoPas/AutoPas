/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "Simulation.h"

#include <autopas/utils/MemoryProfiler.h>
#include <iostream>
#include "PrintableMolecule.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "parsing/MDFlexParser.h"

int main(int argc, char **argv) {
  // start simulation timer
  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> simulation;
  // Parsing
  MDFlexConfig config;

  if (not MDFlexParser::parseInput(argc, argv, config)) {
    exit(-1);
  }
  // make sure sim box is big enough
  config.calcSimulationBox();

  config.print();

  std::cout << std::endl;

  // Initialization
  simulation.initialize(config);
  std::cout << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;

  // Simulation
  std::cout << "Starting simulation... " << std::endl;
  simulation.simulate();
  std::cout << "Simulation done!" << std::endl;

  simulation.printStatistics();

  return EXIT_SUCCESS;
}
