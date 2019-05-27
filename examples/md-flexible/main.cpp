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
void writeVTKFile(string &filename, size_t numParticles, AutoPasTemplate &autopas) {
  std::ofstream vtkFile;
  vtkFile.open(filename);

  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET STRUCTURED_GRID" << endl;
  vtkFile << "DIMENSIONS 1 1 1" << endl;
  vtkFile << "POINTS " << numParticles << " double" << endl;

  for (auto iter = autopas->begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile.close();
}


int main(int argc, char **argv) {
    // Parsing
    MDFlexParser parser;
    if (not parser.parseInput(argc, argv)) {
        exit(-1);
    }
    string logFileName(parser.getLogFileName());
    auto measureFlops(parser.getMeasureFlops());
    auto numIterations(parser.getIterations());
    auto particlesTotal(parser.getParticlesTotal());
    auto particlesPerDim(parser.getParticlesPerDim());
    auto verletRebuildFrequency(parser.getVerletRebuildFrequency());
    auto vtkFilename(parser.getWriteVTK());
    auto logLevel(parser.getLogLevel());

    parser.printConfig();

    //Simulationsdauer ausgerechnet in main.cpp:
    std::chrono::high_resolution_clock::time_point startTotal, stopTotal;
    startTotal = std::chrono::high_resolution_clock::now();

    // select either std::out or a logfile for autopas log output.
    // This does not affect md-flex output.
    std::ofstream logFile;
    std::streambuf *streamBuf;
    if (logFileName.empty()) {
        streamBuf = std::cout.rdbuf();
    } else {
        logFile.open(logFileName);
        streamBuf = logFile.rdbuf();
    }
    std::ostream outputStream(streamBuf);
    long DurationSimulation;


    PrintableMolecule::setEpsilon(1.0);
    PrintableMolecule::setSigma(1.0);
    PrintableMolecule::setMass(1.0);
    //vorübergehend, nicht sicher ob das mit den Werten für OldF Sinn macht
    std::array<double, 3> oldf = {1.0, 1.0, 1.0};
    PrintableMolecule::setOldf(oldf);
    // Initialization
    auto *autopas = new autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> (outputStream);
    autopas::Logger::get()->set_level(logLevel);
    Simulation<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> Simulation(autopas);
    Simulation.initialize(parser);
    //Simulation
    DurationSimulation = Simulation.simulate();
    // @todo -> simulate gibt duration der Simulation zurück

    cout << endl;
    cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
    cout << "sigma  : " << PrintableMolecule::getSigma() << endl << endl;
    if(parser.getGeneratorOption() == MDFlexParser::GeneratorOption::grid){
        particlesTotal=particlesPerDim * particlesPerDim * particlesPerDim;
    }

  if (not vtkFilename.empty()) writeVTKFile(vtkFilename, particlesTotal, autopas);

  // statistics for linked cells
  if (autopas->getContainer()->getContainerType() == autopas::ContainerOption::linkedCells) {
    auto lcContainer = dynamic_cast<autopas::LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> *>(
        autopas->getContainer());
    auto cellsPerDimHalo = lcContainer->getCellBlock().getCellsPerDimensionWithHalo();
    std::array<size_t, 3> cellsPerDim{cellsPerDimHalo[0] - 2, cellsPerDimHalo[1] - 2, cellsPerDimHalo[2] - 2};
    //    auto numCellsHalo = lcContainer->getCells().size();
    auto numCells = cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2];

    cout << "Cells per dimension with Halo: " << cellsPerDimHalo[0] << " x " << cellsPerDimHalo[1] << " x "
         << cellsPerDimHalo[2] << " (Total: " << numCells << ")" << endl;
    cout << "Average Particles per cell: " << (particlesTotal) / (double)numCells << endl;
    cout << endl;
  }

  cout << "Using " << autopas::autopas_get_max_threads() << " Threads" << endl;

  long durationApply = 0;
  unsigned long flopsPerKernelCall = 0;
  cout << "Starting force calculation... " << endl;


   //@todo für statistiken:  flopsPerKernelCall =  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
    //if(MDFlexParser::FunctorOption::lj12_6==parser.getFunctorOption()){
    //  flopsPerKernelCall =
    //          LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
    //} usw für weitere Funktoren


  stopTotal = std::chrono::high_resolution_clock::now();
  cout << "Force calculation done!" << endl;

  //  printMolecules(autopas);

  auto durationTotal = std::chrono::duration_cast<std::chrono::microseconds>(stopTotal - startTotal).count();  //@todo ob durationTotal von main oder von Simulation.h genommen wird
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationApplySec = durationApply * 1e-6;

  //time statistics
  cout << DurationSimulation << endl;
  // Statistics
  cout << fixed << setprecision(2);
  cout << endl << "Measurements:" << endl;
  cout << "Time total   : " << durationTotal << " \u03bcs (" << durationTotalSec << "s)" << endl;
  if (numIterations > 0) {
    cout << "One iteration: " << durationApply / numIterations << " \u03bcs (" << durationApplySec / numIterations
         << "s)" << endl;
  }
  auto mfups = particlesTotal * numIterations / durationApplySec * 1e-6;
  cout << "MFUPs/sec    : " << mfups << endl;

  if (measureFlops) {
    FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(
        autopas->getContainer()->getCutoff());
    autopas->iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * numIterations;
    // approximation for flops of verlet list generation
    if (autopas->getContainer()->getContainerType() == autopas::ContainerOption::verletLists)
      flops +=
          flopCounterFunctor.getDistanceCalculations() *
          FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
          floor(numIterations / verletRebuildFrequency);

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationApplySec << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }

  if (not logFileName.empty()) {
    logFile.close();
  }

  return EXIT_SUCCESS;
}
