/*
 * mdutils.h
 *
 *  Created on: 18 Jan 2018
 *      Author: tchipevn
 */

#ifndef EXAMPLES_MD_MDUTILS_H_
#define EXAMPLES_MD_MDUTILS_H_

#include "autopas.h"
#include <iostream>
#include <cstdlib>

using namespace std;
using namespace autopas;

class PrintableMolecule : public MoleculeLJ {
public:
	PrintableMolecule() : MoleculeLJ() {}
	PrintableMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i) : MoleculeLJ(r, v, i) {}
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

double fRand(double fMin, double fMax) {
	double f = static_cast<double>(rand()) / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

std::array<double, 3> randomPosition(const std::array<double, 3>& boxMin, const std::array<double, 3>& boxMax) {
	std::array<double, 3> r;
	for (int d = 0; d < 3; ++d) {
		r[d] = fRand(boxMin[d], boxMax[d]);
	}
	return r;
}

template<class Molecule, class MoleculeCell>
void fillContainerWithMolecules(int numMolecules, ParticleContainer<Molecule, MoleculeCell>* cont) {
	srand(42); // fixed seedpoint

	std::array<double ,3> boxMin(cont->getBoxMin()), boxMax(cont->getBoxMax());

	for (int i = 0; i < numMolecules; ++i) {
		unsigned long id = static_cast<unsigned long >(i);
		PrintableMolecule m(randomPosition(boxMin, boxMax), {0., 0., 0.}, id);
		cont->addParticle(m);
	}

//	for (auto it = cont->begin(); it.isValid(); ++it) {
//		it->print();
//	}

}


#endif /* EXAMPLES_MD_MDUTILS_H_ */
