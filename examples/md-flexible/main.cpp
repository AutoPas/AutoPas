
#include <AutoPas.h>
#include <chrono>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/GridGenerator.h"
#include "../md/mdutils.h"  // includes autopas.h
#include "MDFlexParser.h"

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
void initContainer(autopas::ContainerOptions containerOption,
                   AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> &autopas, size_t particlesPerDim,
                   double particelSpacing, double cutoff) {
  std::array<double, 3> boxMax({(particlesPerDim + 1.0) * particelSpacing, (particlesPerDim + 1.0) * particelSpacing,
                                (particlesPerDim + 1.0) * particelSpacing});

  // TODO: Extend example to be tunable and also include traversal choices
  autopas.init(boxMax, cutoff, {containerOption});

  PrintableMolecule dummyParticle;
  GridGenerator::fillWithParticles(autopas, {particlesPerDim, particlesPerDim, particlesPerDim}, dummyParticle,
                                   {particelSpacing, particelSpacing, particelSpacing});
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
  auto particleSpacing(parser.getParticlesSpacing());

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal, startCalc, stopCalc;

  startTotal = std::chrono::high_resolution_clock::now();

  // Initialization
  AutoPas<PrintableMolecule, FullParticleCell<PrintableMolecule>> autopas;

  initContainer(containerChoice, autopas, particlesPerDim, particleSpacing, cutoff);

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);
  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma  : " << PrintableMolecule::getSigma() << endl << endl;

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::setGlobals(cutoff, MoleculeLJ::getEpsilon(),
                                                                                MoleculeLJ::getSigma(), 0.0);
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> functor;

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
