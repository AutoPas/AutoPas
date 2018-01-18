/*
 * main.cpp
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */


#include "autopas.h"
#include <iostream>

using namespace std;
using namespace autopas;

class PrintableMolecule : public MoleculeLJ {
public:
	PrintableMolecule() : MoleculeLJ() {}
	PrintableMolecule(std::array<double, 3> r, unsigned long i) : MoleculeLJ(r, i) {}
	void print() {
		cout << "Molecule with position: ";
		for (auto & r : getR()) {
			cout << r << ", " ;
		}
		cout << "and force: ";

		for (auto & f : getF()) {
			cout << f << ", " ;
		}
		cout << "ID: " << getID();
		cout << endl;
	}
};

// how do we recommend to use this thing?
// the user should create a wrapper around our ParticleContainer? That way they can compile whatever they actually need w.r.t. templates once.
// with a bit of documentation, they can

template<class ParticleCell>
void addAFewParticles(ParticleCell& pc) {
	static int i = 0;
	int iEnd = i + 4;
	for ( ; i < iEnd; ++i) {
		std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
		PrintableMolecule m(arr, static_cast<unsigned long>(i));
		pc.addParticle(m);
	}
}

void testForceLJ();

int main(void) {
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
	p1.print();
	p2.print();
	func.AoSFunctor(p2, p1);
	p1.print();
	p2.print();

	testForceLJ();



	cout << "winter is coming" << endl;
	return EXIT_SUCCESS;
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

	for (auto it = container.begin(); it.isValid(); ++it) {
		it->print();
	}

	cout << "done testing iterate pairwise" << endl;
}
