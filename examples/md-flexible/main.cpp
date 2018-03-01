
#include <chrono>
#include <iostream>
#include "../md/mdutils.h"  // includes autopas.h
#include "MDFlexParser.h"

using namespace std;
using namespace autopas;

void fillContainer(
    ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
    *container,
    size_t particlesPerDim) {
  for (unsigned int i = 0; i < particlesPerDim; ++i) {
    for (unsigned int j = 0; j < particlesPerDim; ++j) {
      for (unsigned int k = 0; k < particlesPerDim; ++k) {
        auto p = PrintableMolecule(
            {(double) k + 1, (double) j + 1, (double) i + 1}, {0, 0, 0},
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
    size_t particlesPerDim, double cutoff) {
  std::array<double, 3> boxMin({0., 0., 0.}),
      boxMax({particlesPerDim + 1.0, particlesPerDim + 1.0,
              particlesPerDim + 1.0});

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

  fillContainer(container, particlesPerDim);
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

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal,
      startApply, stopApply;

  startTotal = std::chrono::high_resolution_clock::now();
  ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
      *container = nullptr;
  initContainer(containerChoice, container, particlesPerDim, cutoff);

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);
  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma  : " << PrintableMolecule::getSigma() << endl;

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::setGlobals(
      10.0, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> functorLJ;

  startApply = std::chrono::high_resolution_clock::now();
  apply(*container, functorLJ, dataLayoutChoice);
  stopApply = std::chrono::high_resolution_clock::now();
  stopTotal = std::chrono::high_resolution_clock::now();

  auto durationTotal = std::chrono::duration_cast<std::chrono::microseconds>(
      stopTotal - startTotal)
      .count();
  auto durationApply = std::chrono::duration_cast<std::chrono::microseconds>(
      stopApply - startApply)
      .count();
  auto durationTotalSec =
      std::chrono::duration_cast<std::chrono::seconds>(stopTotal - startTotal)
          .count();
  auto durationApplySec =
      std::chrono::duration_cast<std::chrono::seconds>(stopApply - startApply)
          .count();

  cout << "Time total: " << durationTotal << " \u03bcs (" << durationTotalSec
       << "s)" << endl;
  cout << "Time apply: " << durationApply << " \u03bcs (" << durationApplySec
       << "s)" << endl;

  return EXIT_SUCCESS;
}
