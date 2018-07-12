/**
 * @file main.cpp
 * @date 23.02.2018
 * @author F. Gratl
 */

#include <chrono>
#include <iostream>
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
void printMolecules(AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas) {
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
void initContainer(autopas::ContainerOptions containerOption, std::vector<autopas::TraversalOptions> traversalOptions,
                   AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas, size_t particlesPerDim,
                   double particelSpacing, double cutoff, double verletSkinRadius, int verletRebuildFrequency) {
  std::array<double, 3> boxMax(
      {(particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing, (particlesPerDim)*particelSpacing});

  autopas.init(boxMax, cutoff, verletSkinRadius, verletRebuildFrequency, {containerOption}, traversalOptions);

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particelSpacing, particelSpacing, particelSpacing},
                                   {particelSpacing / 2, particelSpacing / 2, particelSpacing / 2});
}

int main(int argc, char **argv) {
  // Parsing
  MDFlexParser parser;
  if (!parser.parseInput(argc, argv)) {
    exit(-1);
  }

  auto containerChoice(parser.getContainerOption());
  auto dataLayoutChoice(parser.getDataLayoutOption());
  auto particlesPerDim(parser.getParticlesPerDim());
  auto cutoff(parser.getCutoff());
  auto numIterations(parser.getIterations());
  auto particleSpacing(parser.getParticleSpacing());
  auto traversalOptions(parser.getTraversalOptions());
  auto verletRebuildFrequency(parser.getVerletRebuildFrequency());
  auto verletSkinRadius(parser.getVerletSkinRadius());

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal, startCalc, stopCalc;

  startTotal = std::chrono::high_resolution_clock::now();

  // Initialization
  AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> autopas;

  initContainer(containerChoice, traversalOptions, autopas, particlesPerDim, particleSpacing, cutoff, verletSkinRadius,
                verletRebuildFrequency);

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);
  cout << endl;
  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma  : " << PrintableMolecule::getSigma() << endl << endl;

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::setGlobals(cutoff, MoleculeLJ::getEpsilon(),
                                                                                MoleculeLJ::getSigma(), 0.0);
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> functor;

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
    cout << "Particles per cell: " << (particlesPerDim * particlesPerDim * particlesPerDim) / (double)numCells << endl;
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
  auto mmups = particlesPerDim * particlesPerDim * particlesPerDim * numIterations / durationApplySec * 1e-6;

  // Output
  cout << fixed << setprecision(2);
  cout << endl << "Measurements:" << endl;
  cout << "Time total   : " << durationTotal << " \u03bcs (" << durationTotalSec << "s)" << endl;
  cout << "One iteration: " << durationApply / numIterations << " \u03bcs (" << durationApplySec / numIterations << "s)"
       << endl;
  cout << "GFLOPs       : " << flops * 1e-9 << endl;
  cout << "GFLOPs/sec   : " << flops * 1e-9 / durationApplySec << endl;
  cout << "MMUPs/sec    : " << mmups << endl;
  cout << "Hit rate     : " << flopCounterFunctor.getHitRate() << endl;

  return EXIT_SUCCESS;
}
