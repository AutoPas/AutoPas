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

void measure(int which, int numMolecules, int numIterations);

template <class Container>
void measureContainer(Container *cont, int numMolecules, int numIterations);

int main(int argc, char *argv[]) {
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({10., 10., 10.});
  double cutoff = 1.0;

  //	LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> lc;
  //- need to implement addParticle
  //	VerletLists<PrintableMolecule, FullParticleCell<PrintableMolecule>> vl;
  //- need to implement addParticle
  DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> dir(
      boxMin, boxMax, cutoff);

  PrintableMolecule::setEpsilon(1.0);
  PrintableMolecule::setSigma(1.0);

  cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
  cout << "sigma: " << PrintableMolecule::getSigma() << endl;

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
  if (argc == 4) {
    which = atoi(argv[1]);
    numMols = atoi(argv[2]);
    numIts = atoi(argv[3]);

  } else {
    cout << endl
         << "NEEDS THREE ARGUMENTS: <which> <numMolecules> <numIterations>"
         << endl;
    cout << "running: 0(linked cells), 100, 100" << endl << endl;
  }

  measure(which, numMols, numIts);

  cout << "winter is coming" << endl;
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
  cout << "Number of Molecules: " << numMolecules << endl;
  cout << "Number of Force updates: " << numIterations << endl;
  cout << "Elapsed time: " << elapsedTime << endl;
  cout << "MFUPS: " << MFUPS << endl;
  cout << "FLOPs: " << flops << endl;
  cout << "hit rate: " << flopFunctor.getHitRate() << endl;
  cout << "GFLOP/sec:" << flops / elapsedTime * 1e-9 << endl;

  cout << "measuring done" << endl;
}

void measure(int which, int numMolecules, int numIterations) {
  cout << "measuring" << endl;
  std::array<double, 3> boxMin({0., 0., 0.}), boxMax({5., 5., 5.});
  double cutoff = 1.0;

  LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> lcCont(
      boxMin, boxMax, cutoff);
  DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> dirCont(
      boxMin, boxMax, cutoff);

  fillContainerWithMolecules(numMolecules, &lcCont);
  for (auto it = lcCont.begin(); it.isValid(); ++it) {
    dirCont.addParticle(*it);
  }

  if (which == 0) {
    cout << "LINKED CELLS ************************" << endl;
    measureContainer(&lcCont, numMolecules, numIterations);
    cout << "LINKED CELLS DONE *******************" << endl;
  } else if (which == 1) {
    cout << "DIRECT SUM **************************" << endl;
    measureContainer(&dirCont, numMolecules, numIterations);
    cout << "DIRECT SUM DONE *********************" << endl;
  }
}

void testForceLJ() {
  cout << "testing iterate pairwise" << endl;
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

  cout << "done testing iterate pairwise" << endl;
}
