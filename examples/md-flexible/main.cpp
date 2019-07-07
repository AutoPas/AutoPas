/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include "Simulation.h"

#include <autopas/utils/MemoryProfiler.h>
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

/**
 * Writes a VTK file for the current state of the AutoPas object
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
  // starts
  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> simulation;

  // Parsing
  MDFlexParser parser;
  if (not parser.parseInput(argc, argv)) {
    exit(-1);
  }
  auto vtkFilename(parser.getWriteVTK());
  parser.printConfig();
  cout << endl;

  // Initialization
  simulation.initialize(&parser);

<<<<<<< HEAD
  auto autopas = make_shared<autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>>>(outputStream);

  //@todo übernommen vom merge: -> prüfen
  autopas->setTuningStrategyOption(tuningStrategy);
  autopas->setAllowedCellSizeFactors(cellSizeFactors);
  autopas::Logger::get()->set_level(logLevel);

  // setted default anderen boxMax--> sonst Fehler
  autopas->setBoxMax({2., 2., 2.});
  autopas->init();

  autopas::Logger::get()->set_level(logLevel);
  for (auto iter = autopas->begin(); iter.isValid(); ++iter) {
    iter->toString();
  }

  Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> Simulation(autopas);
  Simulation.initialize(&parser);
  // Simulation
  long durationSimulate = Simulation.simulate();
  long durationPosition = Simulation.getDurationX();
  long durationVelocity = Simulation.getDurationV();
  long durationForce = Simulation.getDurationF();
  cout << endl;
  if (parser.getGeneratorOption() == MDFlexParser::GeneratorOption::grid) {
    particlesTotal = particlesPerDim * particlesPerDim * particlesPerDim;
  }

  if (not vtkFilename.empty()) {
    writeVTKFile(vtkFilename, particlesTotal, autopas);
  }
=======
>>>>>>> fe7b66bb907e9a2e9f96c7d01fbae8565c2ec4cb
  cout << "Using " << autopas::autopas_get_max_threads() << " Threads" << endl;

  // Simulation
  cout << "Starting simulation... " << endl;
  simulation.simulate();
  cout << "Simulation done!" << endl;

  simulation.printStatistics();

  return EXIT_SUCCESS;
}
