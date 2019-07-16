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
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "MDFlexParser.h"
#include "PrintableMolecule.h"  // includes autopas.h
#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/LJFunctorAVX.h"

using namespace std;
using namespace autopas;

/**
 * Prints position and forces of all particles in the autopas object.
 * @tparam AutoPasTemplate Template for the templetized autopas type.
 * @param autopas
 */
template <class AutoPasTemplate>
void printMolecules(AutoPasTemplate &autopas) {
  for (auto particleIterator = autopas.begin(); particleIterator.isValid(); ++particleIterator) {
    particleIterator->print();
  }
}

/** Writes a VTK file for the current state of the AutoPas object
 * @tparam AutoPasTemplate Template for the templetized autopas type.
 * @param filename
 * @param numParticles
 * @param autopas
 */
template <class AutoPasTemplate>
void writeVTKFile(string &filename, AutoPasTemplate &autopas) {
  std::ofstream vtkFile;
  vtkFile.open(filename);

  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET STRUCTURED_GRID" << endl;
  vtkFile << "DIMENSIONS 1 1 1" << endl;
  vtkFile << "POINTS " << autopas.getNumberOfParticles() << " double" << endl;

  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile.close();
}

int main(int argc, char **argv) {
  YAML::Node config = YAML::LoadFile("parsingFile.yaml");
  const int particles_per_dimension = config["particles-per-dimension"].as<int>();
  const auto generator = config["generator"].as<std::string>();
  const auto log_level = config["log_level"].as<std::string>();

  std::cout << "particles-per-dimension=  " << particles_per_dimension << endl;
  cout << "generator coice: " << generator << endl;
  cout << "log_level: " << log_level << endl;

  //  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> simulation;
  //  MDFlexParser parser;
  //  if (not parser.parseInput(argc, argv)) {
  //    exit(-1);
  //  }
  //  parser.printConfig();
  //  simulation.initialize(parser);
  //  simulation.printStatistics();
  // frage FABIO, wenn ich hier manuel den destructor von simlation aufrufe; wieso kriege ich 4 invalid reads(autopas
  // container-traversals usw) und 18 invalid free

  return EXIT_SUCCESS;
}
