/*
 * SingleCellIteratorTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#include "autopas.h"
#include "gtest/gtest.h"

using namespace autopas;

template<class ParticleCell>
void addAFewParticles(ParticleCell& pc) {
	static int i = 0;
	int iEnd = i + 4;
	for ( ; i < iEnd; ++i) {
		std::array<double, 3> arr({static_cast<double>(i), static_cast<double>(i), static_cast<double>(i)});
		MoleculeLJ m(arr, {0.,0.,0.}, static_cast<unsigned long>(i));
		pc.addParticle(m);
	}
}


TEST(SingleCellIteratorTest, testFullParticleCell)
{
	FullParticleCell<MoleculeLJ> fpc;
	addAFewParticles(fpc);


	SingleCellIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> it(&fpc);
//	for (; it.isValid(); ++it) {
//		it->print();
//	}
}
