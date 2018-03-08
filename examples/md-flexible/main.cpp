
#include <chrono>
#include <iostream>
#include "../md/mdutils.h"  // includes autopas.h
#include "MDFlexParser.h"

using namespace std;
using namespace autopas;

void fillContainer(
    ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
    *container,
    size_t particlesPerDim, double particelSpacing) {
  for (unsigned int i = 0; i < particlesPerDim; ++i) {
    for (unsigned int j = 0; j < particlesPerDim; ++j) {
      for (unsigned int k = 0; k < particlesPerDim; ++k) {
        auto p = PrintableMolecule(
            {(k + 1) * particelSpacing,
             (j + 1) * particelSpacing,
             (i + 1) * particelSpacing},
            {0, 0, 0},
            i * particlesPerDim * particlesPerDim + j * particlesPerDim + k);
        container->addParticle(p);
      }
    }
  }
}

void printMolecules(
    ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
    *container) {
  for (auto particleIterator = container->begin(); particleIterator.isValid();
       ++particleIterator) {
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
void initContainer(
    MDFlexParser::ContainerOption containerOption,
    ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>> *
    &container,
    size_t particlesPerDim, double particelSpacing, double cutoff) {
  std::array<double, 3> boxMin({0., 0., 0.}),
      boxMax({(particlesPerDim + 1.0) * particelSpacing,
              (particlesPerDim + 1.0) * particelSpacing,
              (particlesPerDim + 1.0) * particelSpacing
             });

  switch (containerOption) {
    case MDFlexParser::directSum: {
      container =
          new DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>>(
              boxMin, boxMax, cutoff);
      break;
    }
    case MDFlexParser::linkedCells: {
      container = new LinkedCells<PrintableMolecule,
                                  FullParticleCell<PrintableMolecule>>(
          boxMin, boxMax, cutoff);
      break;
    }
    default: {
      cout << "Unknown container Option! " << containerOption << endl;
      exit(1);
    }
  }

  fillContainer(container, particlesPerDim, particelSpacing);
}

void apply(
    ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
    &container,
    Functor<PrintableMolecule, FullParticleCell<PrintableMolecule>> &functor,
    MDFlexParser::DataLayoutOption layoutOption) {
  switch (layoutOption) {
    case MDFlexParser::aos: {
      container.iteratePairwiseAoS(&functor);
      break;
    }
    case MDFlexParser::soa: {
      container.iteratePairwiseSoA(&functor);
      break;
    }
  }
}

int main(int argc, char **argv) {
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

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal,
      startApply, stopApply;

  startTotal = std::chrono::high_resolution_clock::now();
  ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
      *container = nullptr;
  initContainer(containerChoice, container, particlesPerDim, particleSpacing, cutoff);

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);
  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma  : " << PrintableMolecule::getSigma() << endl;

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::setGlobals(
      10.0, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> functor;

  startApply = std::chrono::high_resolution_clock::now();
  for (unsigned int i = 0; i < numIterations; ++i) {
    apply(*container, functor, dataLayoutChoice);
  }
  stopApply = std::chrono::high_resolution_clock::now();
  stopTotal = std::chrono::high_resolution_clock::now();

  auto durationTotal = std::chrono::duration_cast<std::chrono::microseconds>(
      stopTotal - startTotal)
      .count();
  auto durationApply = std::chrono::duration_cast<std::chrono::microseconds>(
      stopApply - startApply)
      .count();
  auto durationTotalSec = durationTotal * 1e-6;
  auto durationApplySec = durationApply * 1e-6;

  FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> flopCounterFunctor(container->getCutoff());
  apply(*container, flopCounterFunctor, dataLayoutChoice);
  auto flops = flopCounterFunctor.getFlops(functor.getNumFlopsPerKernelCall()) * numIterations;
  auto mmups = particlesPerDim * particlesPerDim * particlesPerDim * numIterations / durationApplySec * 1e-6;

  cout << fixed << setprecision(2);
  cout << endl << "Measurements:" << endl;
  cout << "Time total: " << durationTotal << " \u03bcs ("
       << durationTotalSec << "s)" << endl;
  cout << "Time apply: " << durationApply / numIterations << " \u03bcs ("
       << durationApplySec / numIterations << "s)" << endl;
  cout << "GFLOPs    : " << flops * 1e-9 << endl;
  cout << "GFLOPs/sec: " << flops * 1e-9 / durationApplySec << endl;
  cout << "MMUPs/sec : " << mmups << endl;
  cout << "Hit rate  : " << flopCounterFunctor.getHitRate() << endl;

  delete container;

  return EXIT_SUCCESS;
}
