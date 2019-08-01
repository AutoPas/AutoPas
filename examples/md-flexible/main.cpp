/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "Simulation.h"

#include <autopas/utils/MemoryProfiler.h>
#include <yaml-cpp/yaml.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "Objects.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "YamlParser.h"
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"

int main(int argc, char **argv) {
  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> simulation;
  auto parser = std::make_shared<YamlParser>();
  if (not parser->parseInput(argc, argv)) {
    exit(-1);
  }
  auto vtkFilename(parser->getWriteVtk());
  parser->printConfig();
  std::cout << std::endl;

  // Initialization
  simulation.initialize(parser);
  std::cout << "Using " << autopas::autopas_get_max_threads() << " Threads" << std::endl;

  // Simulation
  std::cout << "Starting simulation... " << std::endl;
  simulation.simulate();
  std::cout << "Simulation done!" << std::endl;

  simulation.printStatistics();

  return EXIT_SUCCESS;
}
