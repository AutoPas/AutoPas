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
  YamlParser parser;
  //@todo parsing file über die command line übergeben? und default parsingFile definieren
  //@todo catch exception and errors for parser
  std::string filename = "parsingFile.yaml";
  parser.parseInput(filename);
  parser.printConfig();
  simulation.initialize(parser);
  simulation.simulate();
  simulation.printStatistics();

  return EXIT_SUCCESS;
}
