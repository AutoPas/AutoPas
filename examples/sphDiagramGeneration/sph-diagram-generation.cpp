/**
 * @file sph-diagram-generation.cpp
 * @date 25.01.2018
 * @author seckler
 */

#include <autopas/selectors/TraversalSelector.h>
#include <array>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "autopas/autopasIncludes.h"
#include "autopas/sph/autopassph.h"
#include "autopas/utils/Timer.h"

template <class Container, class Functor>
void measureContainer(Container *cont, Functor *func, int numParticles, int numIterations, bool useNewton3);

template <class Container, class Functor, class Traversal>
void measureContainerTraversal(Container *cont, Functor *func, Traversal *traversal, int numParticles,
                               int numIterations);

void addParticles(
    autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> &sph_system,
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
  int containerTypeInt = 0;
  enum ContainerType { linkedCells, directSum, verletLists, verletListsCells } containerType = linkedCells;
  int functorTypeInt = 0;
  enum FunctorType { densityFunctor, hydroForceFunctor } functorType = densityFunctor;

  double skin = 0.;
  int rebuildFrequency = 10;
  bool useNewton3 = true;
  if (argc == 9) {
    numParticles = std::stoi(argv[1]);
    numIterations = std::stoi(argv[2]);
    containerTypeInt = std::stoi(argv[3]);
    functorTypeInt = std::stoi(argv[4]);
    skin = std::stod(argv[5]);
    rebuildFrequency = std::stoi(argv[6]);
    useNewton3 = std::stoi(argv[7]);
    boxMax[0] = boxMax[1] = boxMax[2] = std::stod(argv[8]);
  } else if (argc == 8) {
    numParticles = std::stoi(argv[1]);
    numIterations = std::stoi(argv[2]);
    containerTypeInt = std::stoi(argv[3]);
    functorTypeInt = std::stoi(argv[4]);
    skin = std::stod(argv[5]);
    rebuildFrequency = std::stoi(argv[6]);
    useNewton3 = std::stoi(argv[7]);
  } else if (argc == 7) {
    numParticles = std::stoi(argv[1]);
    numIterations = std::stoi(argv[2]);
    containerTypeInt = std::stoi(argv[3]);
    functorTypeInt = std::stoi(argv[4]);
    skin = std::stod(argv[5]);
    rebuildFrequency = std::stoi(argv[6]);
  } else if (argc == 5) {
    numParticles = std::stoi(argv[1]);
    numIterations = std::stoi(argv[2]);
    containerTypeInt = std::stoi(argv[3]);
    functorTypeInt = std::stoi(argv[4]);
  } else if (argc == 4) {
    numParticles = std::stoi(argv[1]);
    numIterations = std::stoi(argv[2]);
    containerTypeInt = std::stoi(argv[3]);
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

  if (containerTypeInt <= verletListsCells) {
    containerType = static_cast<ContainerType>(containerTypeInt);
  } else {
    std::cerr
        << "Error: wrong containerType " << containerTypeInt << std::endl
        << "containerType should be either 0 (linked-cells), 1 (direct sum), 2 (verlet lists) or 3 (verlet lists cells)"
        << std::endl;
    exit(2);
  }

  if (functorTypeInt <= hydroForceFunctor) {
    functorType = static_cast<FunctorType>(functorTypeInt);
  } else {
    std::cerr << "Error: wrong functorType " << functorTypeInt << std::endl
              << "functorType should be either 0 (density functor) or 1 (hydro force functor)" << std::endl;
    exit(2);
  }

  std::cout << "rebuildFrequency " << rebuildFrequency << " currently unused!" << std::endl;

  autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> lcCont(
      boxMin, boxMax, cutoff, skin * cutoff);
  autopas::DirectSum<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> dirCont(
      boxMin, boxMax, cutoff, skin * cutoff);
  autopas::VerletLists<autopas::sph::SPHParticle> verletCont(boxMin, boxMax, cutoff, skin * cutoff);
  autopas::VerletListsCells<autopas::sph::SPHParticle> verletCellCont(boxMin, boxMax, cutoff,
                                                                      autopas::TraversalOption::c08, skin * cutoff);

  addParticles(lcCont, numParticles);

  for (auto it = lcCont.begin(); it.isValid(); ++it) {
    dirCont.addParticle(*it);
    verletCont.addParticle(*it);
    verletCellCont.addParticle(*it);
  }

  if (containerType == linkedCells) {
    if (functorType == densityFunctor) {
      measureContainer(&lcCont, &densfunc, numParticles, numIterations, useNewton3);
    } else if (functorType == hydroForceFunctor) {
      measureContainer(&lcCont, &hydrofunc, numParticles, numIterations, useNewton3);
    } else {
      std::cout << "wrong functor given" << std::endl;
      exit(2);
    }
  } else if (containerType == directSum) {
    if (functorType == densityFunctor) {
      measureContainer(&dirCont, &densfunc, numParticles, numIterations, useNewton3);
    } else if (functorType == hydroForceFunctor) {
      measureContainer(&dirCont, &hydrofunc, numParticles, numIterations, useNewton3);
    } else {
      std::cout << "wrong functor given" << std::endl;
      exit(2);
    }
  } else if (containerType == verletLists) {
    if (functorType == densityFunctor) {
      measureContainer(&verletCont, &densfunc, numParticles, numIterations, useNewton3);
    } else if (functorType == hydroForceFunctor) {
      measureContainer(&verletCont, &hydrofunc, numParticles, numIterations, useNewton3);
    } else {
      std::cout << "wrong functor given" << std::endl;
      exit(2);
    }
  } else if (containerType == verletListsCells) {
    if (functorType == densityFunctor) {
      measureContainer(&verletCellCont, &densfunc, numParticles, numIterations, useNewton3);
    } else if (functorType == hydroForceFunctor) {
      measureContainer(&verletCellCont, &hydrofunc, numParticles, numIterations, useNewton3);
    } else {
      std::cout << "wrong functor given" << std::endl;
      exit(2);
    }
  } else {
    std::cout << "invalid container id" << std::endl;
    exit(3);
  }
}

template <class Container, class Functor>
void measureContainer(Container *cont, Functor *func, int numParticles, int numIterations, bool useNewton3) {
  using CellType = autopas::FullParticleCell<autopas::sph::SPHParticle>;
  // initialize to dummy values
  auto traversalInfo = cont->getTraversalSelectorInfo();
  std::cout << "Cells: " << traversalInfo.dims[0] << " x " << traversalInfo.dims[1] << " x " << traversalInfo.dims[2]
            << std::endl;
  auto traversalType = autopas::TraversalOption(-1);
  switch (cont->getContainerType()) {
    case autopas::ContainerOption::linkedCells: {
      traversalType = autopas::TraversalOption::c08;
      break;
    }
    case autopas::ContainerOption::directSum: {
      traversalType = autopas::TraversalOption::directSumTraversal;
      //@todo @reviewer ds traversal was previously always with N3 off, ignoring the function arg. Is this intended?
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
    default: {}
  }
  auto traversal = autopas::TraversalSelector<CellType>::template generateTraversal<Functor>(
      traversalType, *func, traversalInfo, autopas::DataLayoutOption::aos,
      useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled);
  measureContainerTraversal(cont, func, dynamic_cast<autopas::CellPairTraversal<CellType> *>(traversal.get()),
                            numParticles, numIterations);
}

template <class Container, class Functor, class Traversal>
void measureContainerTraversal(Container *cont, Functor *func, Traversal *traversal, int numParticles,
                               int numIterations) {
  // autopas::FlopCounterFunctor<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>>
  //    flopFunctor(cont->getCutoff());

  autopas::utils::Timer t;

  // cont->iteratePairwiseAoS(&flopFunctor);
  // double flopsPerIteration = flopFunctor.getFlops(func.getNumFlopsPerKernelCall());

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(traversal);

  double elapsedTime = t.stop();

  // double flops = flopsPerIteration * numIterations;

  double MFUPS_aos = numParticles * numIterations / elapsedTime * 1e-6;

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(traversal);

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
