/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include <iostream>

#include "PrintableMolecule.h"
#include "Simulation.h"
#include "autopas/molecularDynamics/autopasmd.h"
#include "autopas/utils/MemoryProfiler.h"
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

  std::cout << config;

  // Initialization
  simulation.initialize(config);
  std::cout << std::endl << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;

  // Simulation
  std::cout << "Starting simulation... " << std::endl;
  simulation.simulate();
  std::cout << "Simulation done!" << std::endl << std::endl;

  simulation.printStatistics();

  // print config.yaml file of current run
  auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::ostringstream nowStrStr;
  tm unused;
  nowStrStr << std::put_time(localtime_r(&now, &unused), "%Y-%m-%d_%H-%M-%S");
  std::ofstream configFileEnd("MDFlex_end_" + nowStrStr.str() + ".yaml");
  configFileEnd << config;
  configFileEnd.close();

  return EXIT_SUCCESS;
}
