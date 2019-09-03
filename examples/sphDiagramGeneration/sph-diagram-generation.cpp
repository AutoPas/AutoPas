/**
 * @file sph-diagram-generation.cpp
 * @date 25.01.2018
 * @author seckler
 */

#include <autopas/AutoPas.h>
#include <autopas/sph/autopassph.h>
#include <array>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"

template <class Container, class Functor>
void measureContainer(Container *cont, Functor *func, int numParticles, int numIterations);

void addParticles(
    autopas::AutoPas<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> &sph_system,
    int numParticles) {
  // Place SPH particles

  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(sph_system.getBoxMin()), boxMax(sph_system.getBoxMax());

  for (int i = 0; i < numParticles; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::sph::SPHParticle particle(RandomGenerator::randomPosition(boxMin, boxMax), {0., 0., 0.}, id, 0.75, 0.012,
                                       0.);
    // autopas::sph::SPHParticle ith(randomPosition(boxMin, boxMax), {0, 0, 0},
    // i++, 0.75, 0.012, 0. );
    sph_system.addParticle(particle);
  }

  //	for (auto it = cont->begin(); it.isValid(); ++it) {
  //		it->print();
  //	}

  for (auto part = sph_system.begin(); part.isValid(); ++part) {
    part->setMass(part->getMass() * boxMax[0] * boxMax[1] * boxMax[2] / (double)(numParticles));
  }
  // std::cout << "# of ptcls is... " << numParticles << std::endl;

  // Set the end time
}

int main(int argc, char *argv[]) {
  autopas::Logger::create();
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 0.15;
  boxMax[1] = boxMax[2] = boxMax[0] / 1.0;
  double cutoff = .03;

  autopas::sph::SPHCalcDensityFunctor densfunc;
  autopas::sph::SPHCalcHydroForceFunctor hydrofunc;

  int numParticles;
  int numIterations;
  std::set<autopas::ContainerOption> containerOptions;
  int functorTypeInt = 0;
  enum FunctorType { densityFunctor, hydroForceFunctor } functorType = densityFunctor;

  double skin = 0.;
  int rebuildFrequency = 10;
  bool useNewton3 = true;
  if (argc >= 9) {
    boxMax[0] = boxMax[1] = boxMax[2] = std::stod(argv[8]);
  }
  if (argc >= 8) {
    useNewton3 = std::stoi(argv[7]);
  }
  if (argc >= 7) {
    rebuildFrequency = std::stoi(argv[6]);
    skin = std::stod(argv[5]);
  }
  if (argc >= 5) {
    functorTypeInt = std::stoi(argv[4]);
  }
  if (argc >= 4) {
    containerOptions = autopas::utils::StringUtils::parseContainerOptions(argv[3]);
    numIterations = std::stoi(argv[2]);
    numParticles = std::stoi(argv[1]);
  } else {
    std::cerr
        << "ERROR: wrong number of arguments given. " << std::endl
        << "sph-diagram-generation requires the following arguments:" << std::endl
        << "numParticles numIterations containerType [functorType [skin rebuildFrequency [useNewton3 [boxSize]]]]:"
        << std::endl
        << std::endl
        << "containerType should be either 0 (linked-cells), 1 (direct sum), 2 (verlet lists) or 3 (verlet lists cells)"
        << std::endl
        << "functorType should be either 0 (density functor) or 1 (hydro force functor)" << std::endl;
    exit(1);
  }

  if (functorTypeInt <= hydroForceFunctor) {
    functorType = static_cast<FunctorType>(functorTypeInt);
  } else {
    std::cerr << "Error: wrong functorType " << functorTypeInt << std::endl
              << "functorType should be either 0 (density functor) or 1 (hydro force functor)" << std::endl;
    exit(2);
  }

  autopas::AutoPas<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> autoPas;

  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkin(skin * cutoff);
  autoPas.setAllowedContainers(containerOptions);
  autoPas.setAllowedNewton3Options({useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled});
  autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});  // currently aos only!

  auto traversalType = autopas::TraversalOption(-1);
  switch (*containerOptions.begin()) {
    case autopas::ContainerOption::linkedCells: {
      traversalType = autopas::TraversalOption::c08;
      break;
    }
    case autopas::ContainerOption::directSum: {
      traversalType = autopas::TraversalOption::directSumTraversal;
      break;
    }
    case autopas::ContainerOption::verletListsCells: {
      if (useNewton3) {
        traversalType = autopas::TraversalOption::c18;
      } else {
        traversalType = autopas::TraversalOption::c01;
      }
      break;
    }
    case autopas::ContainerOption::verletLists: {
      traversalType = autopas::TraversalOption::verletTraversal;
      break;
    }
    default:
      std::cerr << "Error: containerType " << autopas::utils::StringUtils::to_string(*containerOptions.begin())
                << " not yet supported." << std::endl;
      exit(2);
  }
  autoPas.setAllowedTraversals({traversalType});

  autoPas.init();

  addParticles(autoPas, numParticles);

  if (functorType == densityFunctor) {
    measureContainer(&autoPas, &densfunc, numParticles, numIterations);
  } else if (functorType == hydroForceFunctor) {
    measureContainer(&autoPas, &hydrofunc, numParticles, numIterations);
  } else {
    std::cout << "wrong functor given" << std::endl;
    exit(2);
  }
}

template <class Container, class Functor>
void measureContainer(Container *cont, Functor *func, int numParticles, int numIterations) {

  autopas::utils::Timer t;

  // cont->iteratePairwiseAoS(&flopFunctor);
  // double flopsPerIteration = flopFunctor.getFlops(func.getNumFlopsPerKernelCall());

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(func);

  double elapsedTime = t.stop();

  // double flops = flopsPerIteration * numIterations;

  double MFUPS_aos = numParticles * numIterations / elapsedTime * 1e-6;

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(func);

  elapsedTime = t.stop();

  // double flops = flopsPerIteration * numIterations;

  double MFUPS_soa = numParticles * numIterations / elapsedTime * 1e-6;

  std::cout << numParticles << "\t" << numIterations << "\t" << MFUPS_aos << "\t" << MFUPS_soa;
  // std::cout << numParticles << "\t" << numIterations << "\t" << elapsedTime / numIterations << "\t" << MFUPS;
  // std::cout << "\t" << flops;
  // std::cout << "\t" << flopFunctor.getHitRate();
  // std::cout << "\t" << flops / elapsedTime * 1e-9 << std::endl;

  // std::cout << "measuring done" << std::endl;

  std::cout << std::endl;
}
