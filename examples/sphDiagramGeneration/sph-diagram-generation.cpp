/**
 * @file sph-diagram-generation.cpp
 * @date 25.01.2018
 * @author seckler
 */

#include <array>
#include <iostream>

#include "autopas/AutoPas.h"
#include "autopas/sph/autopassph.h"
#include "autopasTools/generators/RandomGenerator.h"

using Particle = autopas::sph::SPHParticle;
using AutoPasContainer = autopas::AutoPas<Particle>;

template <class Container, class Functor>
void measureContainer(Container *cont, Functor *func, int numParticles, int numIterations);

void addParticles(AutoPasContainer &sph_system, int numParticles) {
  // Place SPH particles

  srand(42);  // fixed seedpoint

  std::array<double, 3> boxMin(sph_system.getBoxMin()), boxMax(sph_system.getBoxMax());

  for (int i = 0; i < numParticles; ++i) {
    auto id = static_cast<unsigned long>(i);
    Particle particle(autopasTools::generators::RandomGenerator::randomPosition(boxMin, boxMax), {0., 0., 0.}, id, 0.75,
                      0.012, 0.);
    sph_system.addParticle(particle);
  }

  for (auto part = sph_system.begin(); part.isValid(); ++part) {
    part->setMass(part->getMass() * boxMax[0] * boxMax[1] * boxMax[2] / (double)(numParticles));
  }

  // Set the end time
}

int main(int argc, char *argv[]) {
  autopas::Logger::create();
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax{};
  boxMax[0] = 0.15;
  boxMax[1] = boxMax[2] = boxMax[0] / 1.0;
  double cutoff = .03;

  autopas::sph::SPHCalcDensityFunctor<Particle> densfunc;
  autopas::sph::SPHCalcHydroForceFunctor<Particle> hydrofunc;

  int numParticles;
  int numIterations;
  std::set<autopas::ContainerOption> containerOptions;
  int functorTypeInt = 0;
  enum FunctorType { densityFunctor, hydroForceFunctor } functorType = densityFunctor;

  double skinPerTimestep = 0.;
  unsigned int rebuildFrequency = 10;
  bool useNewton3 = true;
  try {
    if (argc >= 9) {
      boxMax[0] = boxMax[1] = boxMax[2] = std::stod(argv[8]);
    }
    if (argc >= 8) {
      useNewton3 = std::stoi(argv[7]);
    }
    if (argc >= 7) {
      rebuildFrequency = std::stoi(argv[6]);
      skinPerTimestep = std::stod(argv[5]);
    }
    if (argc >= 5) {
      functorTypeInt = std::stoi(argv[4]);
    }
    if (argc >= 4) {
      containerOptions = autopas::ContainerOption::parseOptions(argv[3]);
      numIterations = std::stoi(argv[2]);
      numParticles = std::stoi(argv[1]);
      if (containerOptions.size() != 1) {
        throw std::runtime_error(
            "Currently choosing of multiple containers is not allowed, please select only one container!");
      }
    } else {
      throw std::runtime_error("too few arguments");
    }
  } catch (const std::exception &e) {
    std::cerr << "ERROR parsing the input arguments: " << e.what() << std::endl
              << "sph-diagram-generation requires the following arguments:" << std::endl
              << "numParticles numIterations containerType [functorType [skinPerTimestep rebuildFrequency [useNewton3 "
                 "[boxSize]]]]:"
              << std::endl
              << std::endl
              << "containerType should be either linked-cells, direct sum, verlet lists or verlet lists cells"
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

  AutoPasContainer autoPas;

  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkinPerTimestep(skinPerTimestep * cutoff / rebuildFrequency);
  autoPas.setVerletRebuildFrequency(rebuildFrequency);
  autoPas.setAllowedContainers(containerOptions);
  autoPas.setAllowedNewton3Options({useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled});
  autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});  // currently aos only!

  autopas::TraversalOption traversalType;
  switch (*containerOptions.begin()) {
    case autopas::ContainerOption::linkedCells: {
      traversalType = autopas::TraversalOption::lc_c08;
      break;
    }
    case autopas::ContainerOption::directSum: {
      traversalType = autopas::TraversalOption::ds_sequential;
      break;
    }
    case autopas::ContainerOption::verletListsCells: {
      if (useNewton3) {
        traversalType = autopas::TraversalOption::lc_c18;
      } else {
        traversalType = autopas::TraversalOption::lc_c01;
      }
      break;
    }
    case autopas::ContainerOption::verletLists: {
      traversalType = autopas::TraversalOption::vl_list_iteration;
      break;
    }
    default:
      std::cerr << "Error: containerType " << containerOptions.begin()->to_string() << " not yet supported."
                << std::endl;
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

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(func);

  double elapsedTime = t.stop();

  double MFUPS_aos = numParticles * numIterations / elapsedTime * 1e-9;

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(func);

  elapsedTime = t.stop();

  double MFUPS_soa = numParticles * numIterations / elapsedTime * 1e-9;

  std::cout << numParticles << "\t" << numIterations << "\t" << MFUPS_aos << "\t" << MFUPS_soa;
  std::cout << std::endl;
}
