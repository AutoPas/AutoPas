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

int main(void) {
//	LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>> lc; - need to implement addParticle
//	VerletLists<PrintableMolecule, FullParticleCell<PrintableMolecule>> vl; - need to implement addParticle
	Direct<PrintableMolecule, FullParticleCell<PrintableMolecule>> dir;

	MoleculeLJ::setEpsilon(1.0 / 7.0);
	MoleculeLJ::setSigma(1.0 / 13.0);

	cout << "epsilon: " << MoleculeLJ::getEpsilon() << endl;
	cout << "sigma: " << MoleculeLJ::getSigma() << endl;

	LJFunctor func;
	func.setGlobals(10.0, 1.0, 1.0, 0.0);
	PrintableMolecule p1({0.0, 0.0, 0.0}, 0);
	PrintableMolecule p2({1.0, 0.0, 0.0}, 1);
	func.AoSFunctor(p1, p2);
	p1.print();
	p2.print();
	func.AoSFunctor(p2, p1);
	p1.print();
	p2.print();



	cout << "winter is coming" << endl;
	return EXIT_SUCCESS;
}
