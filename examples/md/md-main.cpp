/*
 * main.cpp
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#include "autopasIncludes.h"
#include "mdutils.h"
#include "utils/Timer.h"

#include <cstdlib>
#include <iostream>

using namespace std;
using namespace autopas;

void testForceLJ();

void measure(int which, int numMolecules, int numIterations,
             int rebuildFrequency, double skinRadiusToCutoffRatio);

template <class Container>
void measureContainer(Container *cont, int numMolecules, int numIterations);

int main(int argc, char *argv[]) {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 1.0;

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);

//  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
//  cout << "sigma: " << PrintableMolecule::getSigma() << endl;

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>::setGlobals(
      10.0, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);
  PrintableMolecule p1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0);
  PrintableMolecule p2({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1);
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> func;
  func.AoSFunctor(p1, p2);
  //	p1.print();
  //	p2.print();
  func.AoSFunctor(p2, p1);
  //	p1.print();
  //	p2.print();

  testForceLJ();

  int numMols = 100;
  int numIts = 100;
  int which = 0;
  int rebuildFrequency = 1;
  double skinRadius = 0.;
  if (argc == 4) {
    which = atoi(argv[1]);
    numMols = atoi(argv[2]);
    numIts = atoi(argv[3]);

  } else if (argc == 6) {
    which = atoi(argv[1]);
    numMols = atoi(argv[2]);
    numIts = atoi(argv[3]);
    rebuildFrequency = atoi(argv[4]);
  } else {
    cout << endl
         << "NEEDS THREE OR FIVE ARGUMENTS: <which> <numMolecules> "
            "<numIterations> [<rebuildFrequency> <skinradiusToCutoffRatio>]"
         << endl;
    cout << "running: 0(linked cells), 100, 100, 1, 0." << endl << endl;
  }

  measure(which, numMols, numIts, rebuildFrequency, skinRadius);

  return EXIT_SUCCESS;
}

template <class Container>
void measureContainer(Container *cont, int numMolecules, int numIterations) {
  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> func;
  FlopCounterFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>>
      flopFunctor(cont->getCutoff());

  utils::Timer t;

  cont->iteratePairwiseAoS2(&flopFunctor);
  double flopsPerIteration =
      flopFunctor.getFlops(func.getNumFlopsPerKernelCall());

  t.start();
  for (int i = 0; i < numIterations; ++i) {
    cont->iteratePairwiseAoS2(&func);
  }
  double elapsedTime = t.stop();

  double flops = flopsPerIteration * numIterations;

  double MFUPS = numMolecules * numIterations / elapsedTime * 1e-6;
  cout << /*"Number of Molecules: " <<*/ numMolecules<< "\t";// << endl;
  cout << /*"Number of Force updates: " <<*/ numIterations<< "\t";// << endl;
  cout << /*"Elapsed time: " <<*/ elapsedTime<< "\t";// << endl;
  cout << /*"MFUPS: " <<*/ MFUPS<< "\t";// << endl;
  cout << /*"FLOPs: " <<*/ flops<< "\t";// << endl;
  cout << /*"hit rate: " <<*/ flopFunctor.getHitRate()<< "\t";// << endl;
  cout << /*"GFLOP/sec:" <<*/ flops / elapsedTime * 1e-9<< "\t";// << endl;

  cout << endl;
}

void measure(int which, int numMolecules, int numIterations,
             int rebuildFrequency, double skinRadiusToCutoffRatio) {
//  cout << "measuring" << endl;
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({5., 5., 5.});
  double cutoff = 1.0;

  LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> lcCont(
      boxMin, boxMax, cutoff);
  DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> dirCont(
      boxMin, boxMax, cutoff);
  VerletLists<PrintableMolecule, FullParticleCell<PrintableMolecule>>
      verletListsCont(boxMin, boxMax, cutoff, cutoff * skinRadiusToCutoffRatio,
                      rebuildFrequency);

  fillContainerWithMolecules(numMolecules, &lcCont);
  for (auto it = lcCont.begin(); it.isValid(); ++it) {
    dirCont.addParticle(*it);
    verletListsCont.addParticle(*it);
  }

  if (which == 0) {
//    cout << "LINKED CELLS ************************" << endl;
    measureContainer(&lcCont, numMolecules, numIterations);
//    cout << "LINKED CELLS DONE *******************" << endl;
  } else if (which == 1) {
//    cout << "DIRECT SUM **************************" << endl;
    measureContainer(&dirCont, numMolecules, numIterations);
//    cout << "DIRECT SUM DONE *********************" << endl;
  } else if (which == 2) {
//    cout << "VERLET LISTS **************************" << endl;
    measureContainer(&verletListsCont, numMolecules, numIterations);
//    cout << "VERLET LISTS DONE *********************" << endl;
  }
}

void testForceLJ() {
  //cout << "testing iterate pairwise" << endl;
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 1.0;

  DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> container(
      boxMin, boxMax, cutoff);
  PrintableMolecule p1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0);
  PrintableMolecule p2({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1);
  PrintableMolecule p3({0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 2);
  PrintableMolecule p4({1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 3);
  container.addParticle(p1);
  container.addParticle(p2);
  container.addParticle(p3);
  container.addParticle(p4);

  LJFunctor<PrintableMolecule, FullParticleCell<PrintableMolecule>> func;
  container.iteratePairwiseAoS2(&func);

  //	for (auto it = container.begin(); it.isValid(); ++it) {
  //		it->print();
  //	}

  //cout << "done testing iterate pairwise" << endl;
}
