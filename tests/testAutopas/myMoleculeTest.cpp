/*
 * myMoleculeTest.cpp
 *
 *  Created on: Jan 18, 2018
 *      Author: seckler
 */
#include "autopas.h"
#include "gtest/gtest.h"

#include <iostream>

using namespace autopas;
class MyMolecule : public Particle {
public:
	MyMolecule() : Particle(), _myvar(0) {}
	MyMolecule(std::array<double, 3> r, std::array<double, 3> v, unsigned long i, int myvar) : Particle(r, v, i), _myvar(myvar) {}
	void print() {
		std::cout << "Molecule with position: ";
		for (auto & r : getR()) {
			std::cout << r << ", " ;
		}
		std::cout << "and force: ";

		for (auto & f : getF()) {
			std::cout << f << ", " ;
		}
		std::cout << "ID: " << getID();
		std::cout << " myvar: " << _myvar << std::endl;
	}

private:
	int _myvar;
};

template<class ParticleCell>
void addAFewParticles(ParticleCell& pc) {
	static int i = 0;
	int iEnd = i + 4;
	for ( ; i < iEnd; ++i) {
		std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
		MyMolecule m(arr, {0.,0.,0.}, static_cast<unsigned long>(i), i);
		pc.addParticle(m);
	}
}


TEST(myMoleculeTest, testFullParticleCell)
{
	FullParticleCell<MyMolecule> fpc;
	addAFewParticles(fpc);


	SingleCellIterator<MyMolecule, FullParticleCell<MyMolecule>> it(&fpc);
	for (; it.isValid(); ++it) {
		it->print();
	}
}
