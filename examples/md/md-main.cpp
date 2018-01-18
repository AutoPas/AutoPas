/*
 * main.cpp
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */


#include "autopas.h"
#include "mdutils.h"
#include "utils/Timer.h"

#include <iostream>
#include <cstdlib>

using namespace std;
using namespace autopas;

void testForceLJ();

void measureDirect(int numMolecules, int numIterations);

int main(int argc, char * argv[]) {
//	LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> lc; - need to implement addParticle
//	VerletLists<PrintableMolecule, FullParticleCell<PrintableMolecule>> vl; - need to implement addParticle
	DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> dir;

	PrintableMolecule::setEpsilon(1.0);
	PrintableMolecule::setSigma(1.0);

	cout << "epsilon: " << PrintableMolecule::getEpsilon() << endl;
	cout << "sigma: " << PrintableMolecule::getSigma() << endl;

	LJFunctor<PrintableMolecule>::setGlobals(10.0, MoleculeLJ::getEpsilon(), MoleculeLJ::getSigma(), 0.0);
	PrintableMolecule p1({0.0, 0.0, 0.0}, 0);
	PrintableMolecule p2({1.0, 0.0, 0.0}, 1);
	LJFunctor<PrintableMolecule> func;
	func.AoSFunctor(p1, p2);
//	p1.print();
//	p2.print();
	func.AoSFunctor(p2, p1);
//	p1.print();
//	p2.print();

	testForceLJ();

	int numMols = 100;
	int numIts = 100;
	if (argc == 3) {
		numMols = atoi(argv[1]);
		numIts = atoi(argv[2]);
	}

	measureDirect(numMols, numIts);

	cout << "winter is coming" << endl;
	return EXIT_SUCCESS;
}

void measureDirect(int numMolecules, int numIterations) {
	cout << "measuring" << endl;
	DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> cont;
	cont.init();
	std::array<double, 3> boxMin({0.0, 0.0, 0.0});
	std::array<double, 3> boxMax({10.0, 10.0, 10.0});
	cont.setBoxMin(boxMin);
	cont.setBoxMax(boxMax);
	fillContainerWithMolecules(numMolecules, &cont);

	LJFunctor<PrintableMolecule> func;

	utils::Timer t;

	t.start();
	for (int i = 0; i < numIterations; ++i) {
		cont.iteratePairwise(&func);
	}
	double elapsedTime = t.stop();

	double MFUPS = numMolecules * numIterations / elapsedTime * 1e-6;
	cout << "Number of Molecules: " << numMolecules << endl;
	cout << "Number of Force updates: " << numIterations << endl;
	cout << "Elapsed time: " << elapsedTime << endl;
	cout << "MFUPS: " << MFUPS << endl;

//	for (auto it = cont.begin(); it.isValid(); ++it) {
//		it->print();
//	}
	// print one molecule, to be sure compiler isnt optimising stuff away

	int i = 0;
	for (auto it = cont.begin(); it.isValid(); ++it) {
		if (i == 2) {
			break;
		}
		it->print();
		++i;
	}

	cout << "measuring done" << endl;
}

void testForceLJ() {
	cout << "testing iterate pairwise" << endl;

	DirectSum<PrintableMolecule, FullParticleCell<PrintableMolecule>> container;
	container.init();
	PrintableMolecule p1({0.0, 0.0, 0.0}, 0);
	PrintableMolecule p2({1.0, 0.0, 0.0}, 1);
	PrintableMolecule p3({0.0, 1.0, 0.0}, 2);
	PrintableMolecule p4({1.0, 1.0, 0.0}, 3);
	container.addParticle(p1);
	container.addParticle(p2);
	container.addParticle(p3);
	container.addParticle(p4);

	LJFunctor<PrintableMolecule> func;
	container.iteratePairwise(&func);

//	for (auto it = container.begin(); it.isValid(); ++it) {
//		it->print();
//	}

	cout << "done testing iterate pairwise" << endl;
}
