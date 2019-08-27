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
 * @brief Constructs a container and fills it with particles.
 *
 * According to the options passed, a %DirectSum or %'LinkedCells' container is
 * built. It consists of %`FullParticleCells` and is filled with
 * `PrintableMolecules`. The particles are aligned on a cuboid grid.
 *
 * @param autopas AutoPas object that should be initialized
 * @param particlesPerDim Number of desired particles per dimension.
 * @param particleSpacing Space between two particles along each axis of space.
 */
void initContainerGrid(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                       size_t particlesPerDim, double particleSpacing) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax(
      {(particlesPerDim)*particleSpacing, (particlesPerDim)*particleSpacing, (particlesPerDim)*particleSpacing});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particleSpacing, particleSpacing, particleSpacing},
                                   {particleSpacing / 2, particleSpacing / 2, particleSpacing / 2});
}

void initContainerGauss(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                        double boxLength, size_t numParticles, double distributionMean, double distributionStdDev) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  GaussianGenerator::fillWithParticles(autopas, numParticles, dummyParticle, distributionMean, distributionStdDev);
}

void initContainerUniform(autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas,
                          double boxLength, size_t numParticles) {
  std::array<double, 3> boxMin({0., 0., 0.});
  std::array<double, 3> boxMax({boxLength, boxLength, boxLength});

  autopas.setBoxMin(boxMin);
  autopas.setBoxMax(boxMax);

  autopas.init();

  PrintableMolecule dummyParticle;
  RandomGenerator::fillWithParticles(autopas, dummyParticle, numParticles);
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

  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    vtkFile << pos[0] << " " << pos[1] << " " << pos[2] << endl;
  }

  vtkFile.close();
}

/**
 * This function is needed to create functors with the actual type through templates.
 * @tparam FunctorChoice
 * @tparam AutoPasTemplate
 * @param autopas
 * @param cutoff
 * @param numIterations
 * @return Time for all calculation iterations in microseconds.
 */
template <class FunctorChoice, class AutoPasTemplate>
long calculate(AutoPasTemplate &autopas, double cutoff, size_t numIterations) {
  auto functor = FunctorChoice(cutoff, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);

  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;

  startCalc = std::chrono::high_resolution_clock::now();

  // actual Calculation
  for (unsigned int i = 0; i < numIterations; ++i) {
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      cout << "Iteration " << i << endl;
      cout << "Current Memory usage: " << autopas::memoryProfiler::currentMemoryUsage() << " kB" << endl;
    }
    autopas.iteratePairwise(&functor);
  }
  stopCalc = std::chrono::high_resolution_clock::now();

  // TODO: SPEED MULTIPLIER for performance testing purposes ONLY
  float mult = 1.5;
  auto durationCalc = mult * std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
  return durationCalc;
}

int main(int argc, char **argv) {
  // Parsing
  MDFlexParser parser;
  if (not parser.parseInput(argc, argv)) {
    exit(-1);
  }

  auto boxLength(parser.getBoxLength());
  auto containerChoice(parser.getContainerOptions());
  auto selectorStrategy(parser.getSelectorStrategy());
  auto cutoff(parser.getCutoff());
  auto &cellSizeFactors(parser.getCellSizeFactors());
  auto dataLayoutOptions(parser.getDataLayoutOptions());
  auto distributionMean(parser.getDistributionMean());
  auto distributionStdDev(parser.getDistributionStdDev());
  auto functorChoice(parser.getFunctorOption());
  auto generatorChoice(parser.getGeneratorOption());
  auto logLevel(parser.getLogLevel());
  string logFileName(parser.getLogFileName());
  auto measureFlops(parser.getMeasureFlops());
  auto newton3Options(parser.getNewton3Options());
  auto numIterations(parser.getIterations());
  auto particleSpacing(parser.getParticleSpacing());
  auto particlesPerDim(parser.getParticlesPerDim());
  auto particlesTotal(parser.getParticlesTotal());
  auto traversalOptions(parser.getTraversalOptions());
  auto tuningInterval(parser.getTuningInterval());
  auto tuningSamples(parser.getTuningSamples());
  auto tuningStrategy(parser.getTuningStrategyOption());
  auto tuningMaxEvidence(parser.getTuningMaxEvidence());
  auto verletRebuildFrequency(parser.getVerletRebuildFrequency());
  auto verletSkinRadius(parser.getVerletSkinRadius());
  auto vtkFilename(parser.getWriteVTK());

  parser.printConfig();

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
  // Initialization
  autopas::AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> autopas(outputStream);
  autopas::Logger::get()->set_level(logLevel);

  autopas.setCutoff(cutoff);
  autopas.setVerletSkin(verletSkinRadius);
  autopas.setVerletRebuildFrequency(verletRebuildFrequency);
  autopas.setTuningInterval(tuningInterval);
  autopas.setTuningStrategyOption(tuningStrategy);
  autopas.setNumSamples(tuningSamples);
  autopas.setMaxEvidence(tuningMaxEvidence);
  autopas.setSelectorStrategy(selectorStrategy);
  autopas.setAllowedContainers(containerChoice);
  autopas.setAllowedTraversals(traversalOptions);
  autopas.setAllowedDataLayouts(dataLayoutOptions);
  autopas.setAllowedNewton3Options(newton3Options);
  autopas.setAllowedCellSizeFactors(cellSizeFactors);

  switch (generatorChoice) {
    case MDFlexParser::GeneratorOption::grid: {
      initContainerGrid(autopas, particlesPerDim, particleSpacing);
      particlesTotal = particlesPerDim * particlesPerDim * particlesPerDim;
      break;
    }
    case MDFlexParser::GeneratorOption::uniform: {
      initContainerUniform(autopas, boxLength, particlesTotal);
      break;
    }
    case MDFlexParser::GeneratorOption::gaussian: {
      initContainerGauss(autopas, boxLength, particlesTotal, distributionMean, distributionStdDev);
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

  if (not vtkFilename.empty()) {
    writeVTKFile(vtkFilename, particlesTotal, autopas);
  }

  cout << "Using " << autopas::autopas_get_max_threads() << " Threads" << endl;

  long durationApply = 0;
  unsigned long flopsPerKernelCall = 0;
  cout << "Starting force calculation... " << endl;

  switch (functorChoice) {
    case MDFlexParser::FunctorOption::lj12_6: {
      durationApply =
          calculate<LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>>(autopas, cutoff, numIterations);
      flopsPerKernelCall =
          LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
    case MDFlexParser::FunctorOption::lj12_6_AVX: {
      durationApply = calculate<LJFunctorAVX<PrintableMolecule, FullParticleCell<PrintableMolecule>>>(autopas, cutoff,
                                                                                                      numIterations);
      flopsPerKernelCall =
          LJFunctorAVX<PrintableMolecule, FullParticleCell<PrintableMolecule>>::getNumFlopsPerKernelCall();
      break;
    }
  }

  stopTotal = std::chrono::high_resolution_clock::now();
  cout << "Force calculation done!" << endl;

  //  printMolecules(autopas);

  auto durationTotal = std::chrono::duration_cast<std::chrono::microseconds>(stopTotal - startTotal).count();
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationApplySec = durationApply * 1e-6;

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
    FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(autopas.getCutoff());
    autopas.iteratePairwise(&flopCounterFunctor);

    auto flops = flopCounterFunctor.getFlops(flopsPerKernelCall) * numIterations;
    // approximation for flops of verlet list generation
    if (autopas.getContainerType() == autopas::ContainerOption::verletLists)
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
