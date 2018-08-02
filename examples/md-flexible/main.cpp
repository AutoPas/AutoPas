/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include <chrono>
#include <fstream>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GaussianGenerator.h"
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../md/mdutils.h"  // includes autopas.h
#include "MDFlexParser.h"
#include "autopas/AutoPas.h"

using namespace std;
using namespace autopas;

/**
 * Prints position and forces of all particels in the autopas object.
 * @param autopas
 */
void printMolecules(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas) {
  for (auto particleIterator = autopas.begin(); particleIterator.isValid(); ++particleIterator) {
    particleIterator->print();
  }
}

/**
 * @brief Constructs a container and fills it with particles.
 *
 * According to the options passed, a %DirectSum or %'LinkedCells' container is
 * built. It consists of %`FullParticleCells` and is filled with
 * `PrintableMolecules`. The particles are aligned on a cuboid grid.
 *
 * @param containerOption Which container type should be built.
 * @param container Pointer to where the container should be built.
 * @param particlesPerDim Number of desired particles per dimension.
 * @param cutoff Cutoff radius to use. Affects number and size of cells for e.g.
 * LinkedCells.
 */
void initContainerGrid(autopas::ContainerOptions containerOption,
                       std::vector<autopas::TraversalOptions> traversalOptions,
                       autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                       size_t particlesPerDim, double particelSpacing, double cutoff, double verletSkinRadius,
                       int verletRebuildFrequency) {
  std::array<double, 3> boxMax(
      {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

  autopas.init(boxMax, cutoff, verletSkinRadius, verletRebuildFrequency, {containerOption}, traversalOptions);

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particelSpacing, particelSpacing, particelSpacing},
                                   {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

void initContainerGauss(autopas::ContainerOptions containerOption,
                        std::vector<autopas::TraversalOptions> traversalOptions,
                        autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                        double boxLength, size_t numParticles, double distributionMean, double distributionStdDev,
                        double cutoff, double verletSkinRadius, int verletRebuildFrequency) {
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.init(boxMax, cutoff, verletSkinRadius, verletRebuildFrequency, {containerOption}, traversalOptions);
  PrintableMolecule dummyParticle;
  GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

void wirteVTKFile(string filename, size_t numParticles,
                  autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas) {
  std::ofstream vtkFile;
  vtkFile.open(filename);

  vtkFile << "# vtk DataFile Version 2.0" << endl;
  vtkFile << "Timestep" << endl;
  vtkFile << "ASCII" << endl;
  vtkFile << "DATASET STRUCTURED_GRID" << endl;
  vtkFile << "DIMENSIONS 1 1 1" << endl;
  vtkFile << "POINTS " << numParticles << " double" << endl;

  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile.close();
}

int main(int argc, char **argv) {
  // Parsing
  MDFlexParser parser;
  if (!parser.parseInput(argc, argv)) {
    exit(-1);
  }

  auto boxLength(parser.getBoxLength());
  auto containerChoice(parser.getContainerOption());
  auto dataLayoutChoice(parser.getDataLayoutOption());
  auto particlesPerDim(parser.getParticlesPerDim());
  auto cutoff(parser.getCutoff());
  auto numIterations(parser.getIterations());
  auto particleSpacing(parser.getParticleSpacing());
  auto traversalOptions(parser.getTraversalOptions());
  auto verletRebuildFrequency(parser.getVerletRebuildFrequency());
  auto verletSkinRadius(parser.getVerletSkinRadius());
  auto generatorChoice(parser.getGeneratorOption());
  auto distributionMean(parser.getDistributionMean());
  auto distributionStdDev(parser.getDistributionStdDev());
  auto vtkFilename(parser.getWriteVTK());

  parser.printConfig();

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal, startCalc, stopCalc;

  startTotal = std::chrono::high_resolution_clock::now();

  // Initialization
  autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> autopas;

  switch (generatorChoice) {
    case MDFlexParser::GeneratorOption::grid: {
      initContainerGrid(containerChoice, traversalOptions, autopas, particlesPerDim, particleSpacing, cutoff,
                        verletSkinRadius, verletRebuildFrequency);
      break;
    }
    case MDFlexParser::GeneratorOption::gaussian: {
      initContainerGauss(containerChoice, traversalOptions, autopas, boxLength,
                         particlesPerDim * particlesPerDim * particlesPerDim, distributionMean, distributionStdDev,
                         cutoff, verletSkinRadius, verletRebuildFrequency);
      break;
    }
    default:
      std::cerr << "Unknown generator choice" << std::endl;
      return -1;
  }

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);
  cout << endl;
  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma  : " << PrintableMolecule::getSigma() << endl << endl;

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::setGlobals(cutoff, MoleculeLJ::getEpsilon(),
                                                                                MoleculeLJ::getSigma(), 0.0);
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> functor;

  if (not vtkFilename.empty()) wirteVTKFile(vtkFilename, particlesPerDim * particlesPerDim * particlesPerDim, autopas);

  // statistics for linked cells
  if (containerChoice == autopas::ContainerOptions::linkedCells) {
    auto lcContainer = dynamic_cast<autopas::LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> *>(
        autopas.getContainer());
    auto cellsPerDimHalo = lcContainer->getCellBlock().getCellsPerDimensionWithHalo();
    std::array<size_t, 3> cellsPerDim{cellsPerDimHalo[0] - 2, cellsPerDimHalo[1] - 2, cellsPerDimHalo[2] - 2};
    auto numCellsHalo = lcContainer->getCells().size();
    auto numCells = cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2];

    cout << "Cells per dimension with Halo: " << cellsPerDimHalo[0] << " x " << cellsPerDimHalo[1] << " x "
         << cellsPerDimHalo[2] << " (Total: " << numCells << ")" << endl;
    cout << "Average Particles per cell: " << (particlesPerDim * particlesPerDim * particlesPerDim) / (double)numCells
         << endl;
    cout << endl;
  }

  cout << "Using " << autopas::autopas_get_max_threads() << " Threads" << endl;
  cout << "Starting force calculation... " << flush;
  startCalc = std::chrono::high_resolution_clock::now();
  // Calculation
  for (unsigned int i = 0; i < numIterations; ++i) {
    autopas.iteratePairwise(&functor, dataLayoutChoice);
  }
  stopCalc = std::chrono::high_resolution_clock::now();
  stopTotal = std::chrono::high_resolution_clock::now();
  cout << "done!" << endl;

  //  printMolecules(autopas);

  // Statistics
  auto durationTotal = std::chrono::duration_cast<std::chrono::microseconds>(stopTotal - startTotal).count();
  auto durationApply = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationApplySec = durationApply * 1e-6;

  FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(
      autopas.getContainer()->getCutoff());
  autopas.iteratePairwise(&flopCounterFunctor, dataLayoutChoice);

  auto flops = flopCounterFunctor.getFlops(functor.getNumFlopsPerKernelCall()) * numIterations;
  // approximation for flops of verlet list generation
  if (containerChoice == autopas::ContainerOptions::verletLists)
    flops +=
        flopCounterFunctor.getDistanceCalculations() *
        FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
        floor(numIterations / verletRebuildFrequency);
  auto mfups = particlesPerDim * particlesPerDim * particlesPerDim * numIterations / durationApplySec * 1e-6;

  // Output
  cout << fixed << setprecision(2);
  cout << endl << "Measurements:" << endl;
  cout << "Time total   : " << durationTotal << " \u03bcs (" << durationTotalSec << "s)" << endl;
  if (numIterations > 0)
    cout << "One iteration: " << durationApply / numIterations << " \u03bcs (" << durationApplySec / numIterations
         << "s)" << endl;
  cout << "GFLOPs       : " << flops * 1e-9 << endl;
  cout << "GFLOPs/sec   : " << flops * 1e-9 / durationApplySec << endl;
  cout << "MFUPs/sec    : " << mfups << endl;
  cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;

  return EXIT_SUCCESS;
}
