/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

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
 * @param containerOptions Which container type should be built.
 * @param container Pointer to where the container should be built.
 * @param particlesPerDim Number of desired particles per dimension.
 * @param cutoff Cutoff radius to use. Affects number and size of cells for e.g.
 * LinkedCells.
 */
void initContainerGrid(std::vector<autopas::ContainerOptions> &containerOptions,
                       autopas::SelectorStrategy containerSelectorStrategy,
                       std::vector<autopas::TraversalOptions> &traversalOptions,
                       autopas::SelectorStrategy traversalSelectorStrategy,
                       autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                       size_t particlesPerDim, double particelSpacing, double cutoff, double verletSkinRadius,
                       unsigned int verletRebuildFrequency, unsigned int tuningInterval, unsigned int tuningSamples) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax(
      {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

  autopas.init(boxMin, boxMax, cutoff, verletSkinRadius, verletRebuildFrequency, containerOptions, traversalOptions,
               containerSelectorStrategy, traversalSelectorStrategy, tuningInterval, tuningSamples);

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particelSpacing, particelSpacing, particelSpacing},
                                   {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

void initContainerGauss(std::vector<autopas::ContainerOptions> &containerOptions,
                        autopas::SelectorStrategy containerSelectorStrategy,
                        std::vector<autopas::TraversalOptions> &traversalOptions,
                        autopas::SelectorStrategy traversalSelectorStrategy,
                        autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                        double boxLength, size_t numParticles, double distributionMean, double distributionStdDev,
                        double cutoff, double verletSkinRadius, unsigned int verletRebuildFrequency,
                        unsigned int tuningInterval, unsigned int tuningSamples) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.init(boxMin, boxMax, cutoff, verletSkinRadius, verletRebuildFrequency, containerOptions, traversalOptions,
               containerSelectorStrategy, traversalSelectorStrategy, tuningInterval, tuningSamples);

  PrintableMolecule dummyParticle;
  GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

void initContainerUniform(std::vector<autopas::ContainerOptions> &containerOptions,
                          autopas::SelectorStrategy containerSelectorStrategy,
                          std::vector<autopas::TraversalOptions> &traversalOptions,
                          autopas::SelectorStrategy traversalSelectorStrategy,
                          autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                          double boxLength, size_t numParticles, double cutoff, double verletSkinRadius,
                          unsigned int verletRebuildFrequency, unsigned int tuningInterval,
                          unsigned int tuningSamples) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.init(boxMin, boxMax, cutoff, verletSkinRadius, verletRebuildFrequency, containerOptions, traversalOptions,
               containerSelectorStrategy, traversalSelectorStrategy, tuningInterval, tuningSamples);

  PrintableMolecule dummyParticle;
  RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
}

void wirteVTKFile(string &filename, size_t numParticles,
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
  if (not parser.parseInput(argc, argv)) {
    exit(-1);
  }

  auto boxLength(parser.getBoxLength());
  auto containerChoice(parser.getContainerOptions());
  auto containerSelectorStrategy(parser.getContainerSelectorStrategy());
  auto cutoff(parser.getCutoff());
  auto dataLayoutChoice(parser.getDataLayoutOption());
  auto distributionMean(parser.getDistributionMean());
  auto distributionStdDev(parser.getDistributionStdDev());
  auto generatorChoice(parser.getGeneratorOption());
  auto logLevel(parser.getLogLevel());
  auto measureFlops(parser.getMeasureFlops());
  auto numIterations(parser.getIterations());
  auto particleSpacing(parser.getParticleSpacing());
  auto particlesTotal(parser.getParticlesTotal());
  auto particlesPerDim(parser.getParticlesPerDim());
  auto traversalOptions(parser.getTraversalOptions());
  auto traversalSelectorStrategy(parser.getTraversalSelectorStrategy());
  auto tuningInterval(parser.getTuningInterval());
  auto tuningSamples(parser.getTuningSamples());
  auto verletRebuildFrequency(parser.getVerletRebuildFrequency());
  auto verletSkinRadius(parser.getVerletSkinRadius());
  auto vtkFilename(parser.getWriteVTK());

  parser.printConfig();

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal, startCalc, stopCalc;

  startTotal = std::chrono::high_resolution_clock::now();

  // Initialization
  autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> autopas;
  autopas::Logger::get()->set_level(logLevel);
  switch (generatorChoice) {
    case MDFlexParser::GeneratorOption::grid: {
      initContainerGrid(containerChoice, containerSelectorStrategy, traversalOptions, traversalSelectorStrategy,
                        autopas, particlesPerDim, particleSpacing, cutoff, verletSkinRadius, verletRebuildFrequency,
                        tuningInterval, tuningSamples);
      particlesTotal = particlesPerDim * particlesPerDim * particlesPerDim;
      break;
    }
    case MDFlexParser::GeneratorOption::uniform: {
      initContainerUniform(containerChoice, containerSelectorStrategy, traversalOptions, traversalSelectorStrategy,
                           autopas, boxLength, particlesTotal, cutoff, verletSkinRadius, verletRebuildFrequency,
                           tuningInterval, tuningSamples);
      break;
    }
    case MDFlexParser::GeneratorOption::gaussian: {
      initContainerGauss(containerChoice, containerSelectorStrategy, traversalOptions, traversalSelectorStrategy,
                         autopas, boxLength, particlesTotal, distributionMean, distributionStdDev, cutoff,
                         verletSkinRadius, verletRebuildFrequency, tuningInterval, tuningSamples);
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

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> functor(cutoff, MoleculeLJ::getEpsilon(),
                                                                            MoleculeLJ::getSigma(), 0.0);

  if (not vtkFilename.empty()) wirteVTKFile(vtkFilename, particlesTotal, autopas);

  // statistics for linked cells
  if (autopas.getContainer()->getContainerType() == autopas::ContainerOptions::linkedCells) {
    auto lcContainer = dynamic_cast<autopas::LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> *>(
        autopas.getContainer());
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
  cout << "Starting force calculation... " << endl;
  startCalc = std::chrono::high_resolution_clock::now();
  // Calculation
  for (unsigned int i = 0; i < numIterations; ++i) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      cout << "Iteration " << i << endl;
      cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
    }
    autopas.iteratePairwise(&functor, dataLayoutChoice);
  }
  stopCalc = std::chrono::high_resolution_clock::now();
  stopTotal = std::chrono::high_resolution_clock::now();
  cout << "Force calculation done!" << endl;

  //  printMolecules(autopas);

  auto durationTotal = std::chrono::duration_cast<std::chrono::microseconds>(stopTotal - startTotal).count();
  auto durationApply = std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationApplySec = durationApply * 1e-6;

  // Statistics
  cout << fixed << setprecision(2);
  cout << endl << "Measurements:" << endl;
  cout << "Time total   : " << durationTotal << " \u03bcs (" << durationTotalSec << "s)" << endl;
  if (numIterations > 0)
    cout << "One iteration: " << durationApply / numIterations << " \u03bcs (" << durationApplySec / numIterations
         << "s)" << endl;

  if (measureFlops) {
    FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(
        autopas.getContainer()->getCutoff());
    autopas.iteratePairwise(&flopCounterFunctor, dataLayoutChoice);

    auto flops = flopCounterFunctor.getFlops(functor.getNumFlopsPerKernelCall()) * numIterations;
    // approximation for flops of verlet list generation
    if (autopas.getContainer()->getContainerType() == autopas::ContainerOptions::verletLists)
      flops +=
          flopCounterFunctor.getDistanceCalculations() *
          FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::numFlopsPerDistanceCalculation *
          floor(numIterations / verletRebuildFrequency);
    auto mfups = particlesTotal * numIterations / durationApplySec * 1e-6;

    cout << "GFLOPs       : " << flops * 1e-9 << endl;
    cout << "GFLOPs/sec   : " << flops * 1e-9 / durationApplySec << endl;
    cout << "MFUPs/sec    : " << mfups << endl;
    cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;
  }

  return EXIT_SUCCESS;
}
