
#include <chrono>
#include <iostream>
#include "../md/mdutils.h"
#include "MDFlexParser.h"
#include "autopas.h"

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

void initContainer(
    MDFlexParser::ContainerOption containerOption,
    ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>> *
    &container,
    size_t particlesPerDim) {
  std::array<double, 3> boxMin({0., 0., 0.}),
      boxMax({particlesPerDim + 1.0, particlesPerDim + 1.0,
              particlesPerDim + 1.0});
  double cutoff = 1.0;

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

void apply(ParticleContainer<PrintableMolecule,
                             FullParticleCell<PrintableMolecule>> &container,
           Functor<PrintableMolecule> &functor,
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
    cerr << "Could not parse arguments!" << endl;
    exit(-1);
  }

  MDFlexParser::ContainerOption containerChoice(parser.getContainerOption());
  MDFlexParser::DataLayoutOption dataLayoutChoice(parser.getDataLayoutOption());
  size_t particlesPerDim(parser.getParticlesPerDim());

  std::chrono::high_resolution_clock::time_point startTotal, stopTotal, startApply, stopApply;

  startTotal = std::chrono::high_resolution_clock::now();
  ParticleContainer<PrintableMolecule, FullParticleCell<PrintableMolecule>>
      *container = nullptr;
  initContainer(containerChoice, container, particlesPerDim);

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);
  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma  : " << PrintableMolecule::getSigma() << endl;

  LJFunctor<PrintableMolecule>::setGlobals(10.0, MoleculeLJ::getEpsilon(),
                                           MoleculeLJ::getSigma(), 0.0);
  LJFunctor<PrintableMolecule> functorLJ;

  startApply = std::chrono::high_resolution_clock::now();
  apply(*container, functorLJ, dataLayoutChoice);
  stopApply = std::chrono::high_resolution_clock::now();
  stopTotal = std::chrono::high_resolution_clock::now();

//  printMolecules(container);

  auto durationTotal =
      std::chrono::duration_cast<std::chrono::microseconds>(stopTotal - startTotal)
          .count();
  auto durationApply =
      std::chrono::duration_cast<std::chrono::microseconds>(stopApply - startApply)
          .count();

  cout << "Time total: " << durationTotal << " \u03bcs" << endl;
  cout << "Time apply: " << durationApply << " \u03bcs" << endl;

  return EXIT_SUCCESS;
}
