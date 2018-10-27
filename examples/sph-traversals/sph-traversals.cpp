/**
 * @file sph-traversals.cpp
 * @date 20.09.18
 * @author nguyen
 */

#include <array>
#include <iostream>
#include "../md/mdutils.h"
#include "autopas/autopasIncludes.h"
#include "autopas/sph/autopassph.h"
#include "autopas/utils/Timer.h"

template <class Container, class Traversal>
void measureContainer(Container *cont, autopas::sph::SPHCalcDensityFunctor *func, Traversal *traversal,
                      int numParticles, int numIterations, bool useNewton3);

void addParticles(
    autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> &sph_system,
    int numParticles) {
  // Place SPH particles

  srand(10032);  // fixed seedpoint

  std::array<double, 3> boxMin(sph_system.getBoxMin()), boxMax(sph_system.getBoxMax());

  for (int i = 0; i < numParticles; ++i) {
    auto id = static_cast<unsigned long>(i);
    autopas::sph::SPHParticle particle(randomPosition(boxMin, boxMax), {0., 0., 0.}, id, 0.75, 0.012, 0.);
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

  autopas::sph::SPHCalcDensityFunctor func;

  int numParticles = 16;
  int numIterations = 100000;
  int containerTypeInt = 0;
  enum ContainerType { linkedCells, verletListsCells, verletCluster } containerType;
  int traversalInt = 0;
  enum TraversalType { c01, c08, c18, sliced } traversalType;

  bool useNewton3 = true;
  double skin = 0.;
  int rebuildFrequency = 10;
  if (argc == 9) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
    boxMax[0] = boxMax[1] = boxMax[2] = atof(argv[4]);
    traversalInt = atoi(argv[5]);
    useNewton3 = atoi(argv[6]);
    skin = atof(argv[7]);
    rebuildFrequency = atof(argv[8]);
  } else if (argc == 7) {
    numParticles = atoi(argv[1]);
    numIterations = atoi(argv[2]);
    containerTypeInt = atoi(argv[3]);
    boxMax[0] = boxMax[1] = boxMax[2] = atof(argv[4]);
    traversalInt = atoi(argv[5]);
    useNewton3 = atoi(argv[6]);
  } else {
    std::cerr << "ERROR: wrong number of arguments given. " << std::endl
              << "sph-traversals requires the following arguments:" << std::endl
              << "numParticles numIterations containerType boxSize traversal useNewton3 [skin rebuildFrequency]:"
              << std::endl
              << std::endl
              << "containerType should be either 0 (linked-cells), 1 (verlet lists cells), 2 (verlet cluster)"
              << std::endl
              << "traversal should be either 0 (c01), 1 (c08), 2 (c18) or 3 (sliced)" << std::endl;
    exit(1);
  }

  if (containerTypeInt <= verletCluster) {
    containerType = static_cast<ContainerType>(containerTypeInt);
  } else {
    std::cerr << "Error: wrong containerType " << containerTypeInt << std::endl
              << "containerType should be either 0 (linked-cells), 1 (verlet lists cells), 2 (verlet cluster)"
              << std::endl;
    exit(2);
  }

  if (traversalInt <= sliced) {
    traversalType = static_cast<TraversalType>(traversalInt);
  } else {
    std::cerr << "Error: wrong traversalType " << traversalInt << std::endl
              << "traversal should be either 0 (c01), 1 (c08), 2 (c18) or 3 (sliced)" << std::endl;
    exit(2);
  }

  autopas::LinkedCells<autopas::sph::SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> lcCont(
      boxMin, boxMax, cutoff);
  autopas::VerletListsCells<autopas::sph::SPHParticle> verletCellContc08(
      boxMin, boxMax, cutoff, autopas::TraversalOptions::c08, skin * cutoff, rebuildFrequency);
  autopas::VerletListsCells<autopas::sph::SPHParticle> verletCellContc18(
      boxMin, boxMax, cutoff, autopas::TraversalOptions::c18, skin * cutoff, rebuildFrequency);
  autopas::VerletClusterLists<autopas::sph::SPHParticle> verletClusterCont(boxMin, boxMax, cutoff, skin * cutoff,
                                                                           rebuildFrequency);

  addParticles(lcCont, numParticles);

  for (auto it = lcCont.begin(); it.isValid(); ++it) {
    verletCellContc08.addParticle(*it);
    verletCellContc18.addParticle(*it);
    verletClusterCont.addParticle(*it);
  }

  if (containerType == linkedCells) {
    auto dims = lcCont.getCellBlock().getCellsPerDimensionWithHalo();
    std::cout << "Cells: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;

    if (traversalType == c01) {
      if (not useNewton3) {
        C01Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false> traversal(
            dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        std::cout << "c01 does not support newton3" << std::endl;
        exit(3);
      }
    } else if (traversalType == c08) {
      if (useNewton3) {
        C08Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, true>
            traversal(dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        C08Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, false>
            traversal(dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      }
    } else if (traversalType == c18) {
      if (useNewton3) {
        C18Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, true>
            traversal(dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        C18Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, false>
            traversal(dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      }
    } else if (traversalType == sliced) {
      if (useNewton3) {
        SlicedTraversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, true>
            traversal(dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        SlicedTraversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, false>
            traversal(dims, &func);
        measureContainer(&lcCont, &func, &traversal, numParticles, numIterations, useNewton3);
      }
    } else {
      std::cout << "invalid traversal id" << std::endl;
      exit(3);
    }
  } else if (containerType == verletListsCells) {
    auto dims = verletCellContc18.getCellsPerDimension();
    std::cout << "Cells: " << dims[0] << " x " << dims[1] << " x " << dims[2] << std::endl;

    if (traversalType == c01) {
      if (not useNewton3) {
        C01Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false> traversal(
            dims, &func);
        measureContainer(&verletCellContc08, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        std::cout << "c01 does not support newton3" << std::endl;
        exit(3);
      }
    } else if (traversalType == c08) {
      std::cout << "c08 not implemented yet" << std::endl;
      exit(3);
    } else if (traversalType == c18) {
      if (useNewton3) {
        C18Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, true>
            traversal(dims, &func);
        measureContainer(&verletCellContc18, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        C18Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, false>
            traversal(dims, &func);
        measureContainer(&verletCellContc18, &func, &traversal, numParticles, numIterations, useNewton3);
      }
    } else if (traversalType == sliced) {
      if (useNewton3) {
        SlicedTraversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, true>
            traversal(dims, &func);
        measureContainer(&verletCellContc08, &func, &traversal, numParticles, numIterations, useNewton3);
      } else {
        SlicedTraversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false, false>
            traversal(dims, &func);
        measureContainer(&verletCellContc08, &func, &traversal, numParticles, numIterations, useNewton3);
      }
    } else {
      std::cout << "invalid traversal id" << std::endl;
      exit(3);
    }
  } else if (containerType == verletCluster) {
    if (traversalType == c01) {
      if (not useNewton3) {
        C01Traversal<FullParticleCell<autopas::sph::SPHParticle>, autopas::sph::SPHCalcDensityFunctor, false>
            dummyTraversal({0, 0, 0}, &func);
        measureContainer(&verletClusterCont, &func, &dummyTraversal, numParticles, numIterations, useNewton3);
      } else {
        std::cout << "c01 does not support newton3" << std::endl;
        exit(3);
      }
    } else {
      std::cout << "traversal invalid or not implemented yet" << std::endl;
      exit(3);
    }
  } else {
    std::cout << "invalid container id" << std::endl;
    exit(4);
  }
}

template <class Container, class Traversal>
void measureContainer(Container *cont, autopas::sph::SPHCalcDensityFunctor *func, Traversal *traversal,
                      int numParticles, int numIterations, bool useNewton3) {
  autopas::utils::Timer t;

  t.start();
  for (int i = 0; i < numIterations; ++i) {
    cont->iteratePairwiseAoS(func, traversal, useNewton3);
  }
  double elapsedTime = t.stop();
  double MFUPS = numParticles * numIterations / elapsedTime * 1e-6;

  std::cout << "MFUPS: " << MFUPS << std::endl;
}
