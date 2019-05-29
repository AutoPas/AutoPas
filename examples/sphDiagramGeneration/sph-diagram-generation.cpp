/**
 * @file sph-diagram-generation.cpp
 * @date 25.01.2018
 * @author seckler
 */

#include <array>
#include <iostream>
#include "../../tests/testAutopas/testingHelpers/RandomGenerator.h"
#include "autopas/autopasIncludes.h"
#include "autopas/containers/cellPairTraversals/DummyTraversal.h"
#include "autopas/sph/autopassph.h"
#include "autopas/utils/Timer.h"

template <class Container, class Functor>
void measureContainer(Container *cont, Functor *func, int numParticles, int numIterations, bool useNewton3);

template <class Container, class Functor, class Traversal>
void measureContainerTraversal(Container *cont, Functor *func, Traversal *traversal, int numParticles,
                               int numIterations, bool useNewton3);

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

  int numParticles = 16;
  int numIterations = 100000;
  int containerTypeInt = 0;
  enum ContainerType { linkedCells, directSum, verletLists, verletListsCells } containerType = linkedCells;
  int functorTypeInt = 0;
  enum FunctorType { densityFunctor, hydroForceFunctor } functorType = densityFunctor;

  double skin = 0.;
  int rebuildFrequency = 10;
  bool useNewton3 = true;
  if (argc == 9) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
    functorTypeInt = atoi(argv[4]);
    skin = atof(argv[5]);
    rebuildFrequency = atof(argv[6]);
    useNewton3 = atoi(argv[7]);
    boxMax[0] = boxMax[1] = boxMax[2] = atof(argv[8]);
  } else if (argc == 8) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
    functorTypeInt = atoi(argv[4]);
    skin = atof(argv[5]);
    rebuildFrequency = atof(argv[6]);
    useNewton3 = atoi(argv[7]);
  } else if (argc == 7) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
    functorTypeInt = atoi(argv[4]);
    skin = atof(argv[5]);
    rebuildFrequency = atof(argv[6]);
  } else if (argc == 5) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
    functorTypeInt = atoi(argv[4]);
  } else if (argc == 4) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
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

  autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> lcCont(
      boxMin, boxMax, cutoff);
  autopas::DirectSum<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> dirCont(
      boxMin, boxMax, cutoff);
  autopas::VerletLists<autopas::sph::SPHParticle> verletCont(boxMin, boxMax, cutoff, skin * cutoff, rebuildFrequency);
  autopas::VerletListsCells<autopas::sph::SPHParticle> verletCellCont(
      boxMin, boxMax, cutoff, autopas::TraversalOption::c08, skin * cutoff, rebuildFrequency);

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
  autopas::CellPairTraversal<autopas::FullParticleCell<autopas::sph::SPHParticle>> *traversal;

  switch (cont->getContainerType()) {
    case autopas::ContainerOption::linkedCells: {
      auto dims =
          dynamic_cast<
              autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> *>(
              cont)
              ->getCellBlock()
              .getCellsPerDimensionWithHalo();
      std::cout << "Cells: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;

      if (useNewton3) {
        traversal = new autopas::C08Traversal<autopas::FullParticleCell<autopas::sph::SPHParticle>, Functor,
                                              autopas::DataLayoutOption::aos, true>(dims, func);
      } else {
        traversal = new autopas::C08Traversal<autopas::FullParticleCell<autopas::sph::SPHParticle>, Functor,
                                              autopas::DataLayoutOption::aos, false>(dims, func);
      }

      break;
    }

    case autopas::ContainerOption::directSum: {
      traversal = new autopas::DirectSumTraversal<autopas::FullParticleCell<autopas::sph::SPHParticle>, Functor,
                                                  autopas::DataLayoutOption::aos, false>(func);
      break;
    }
    case autopas::ContainerOption::verletListsCells: {
      auto dims = dynamic_cast<autopas::VerletListsCells<autopas::sph::SPHParticle> *>(cont)->getCellsPerDimension();
      std::cout << "Cells: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;

      if (useNewton3) {
        traversal = new autopas::C18Traversal<autopas::FullParticleCell<autopas::sph::SPHParticle>, Functor,
                                              autopas::DataLayoutOption::aos, true>(dims, func);
      } else {
        traversal = new autopas::C01Traversal<autopas::FullParticleCell<autopas::sph::SPHParticle>, Functor,
                                              autopas::DataLayoutOption::aos, false>(dims, func);
      }
      break;
    }
    default:
      traversal = new autopas::DummyTraversal<autopas::FullParticleCell<autopas::sph::SPHParticle>>({0, 0, 0});
  }

  measureContainerTraversal(cont, func, traversal, numParticles, numIterations, useNewton3);
}

template <class Container, class Functor, class Traversal>
void measureContainerTraversal(Container *cont, Functor *func, Traversal *traversal, int numParticles,
                               int numIterations, bool useNewton3) {
  // autopas::FlopCounterFunctor<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>>
  //    flopFunctor(cont->getInteractionLength());

  autopas::utils::Timer t;

  // cont->iteratePairwiseAoS(&flopFunctor);
  // double flopsPerIteration = flopFunctor.getFlops(func.getNumFlopsPerKernelCall());

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(func, traversal, useNewton3);

  double elapsedTime = t.stop();

  // double flops = flopsPerIteration * numIterations;

  double MFUPS_aos = numParticles * numIterations / elapsedTime * 1e-6;

  t.start();
  for (int i = 0; i < numIterations; ++i) cont->iteratePairwise(func, traversal, useNewton3);

  elapsedTime = t.stop();

  delete traversal;

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
