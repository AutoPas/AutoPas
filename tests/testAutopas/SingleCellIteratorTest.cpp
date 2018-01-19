/*
 * SingleCellIteratorTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#include "SingleCellIteratorTest.h"
#include "autopas.h"
#include "gtest/gtest.h"

using namespace autopas;

void SingleCellIteratorTest::SetUp() {
	for (int i = 0; i < 4; ++i) {
		std::array<double, 3> arr;
		for (auto & a : arr) {
			a = static_cast<double>(i);
		}
		MoleculeLJ m(arr, { 0., 0., 0. }, static_cast<unsigned long>(i));
		_vecOfMolecules.push_back(m);
	}
}

void SingleCellIteratorTest::fillWithParticles(
		autopas::ParticleCell<autopas::MoleculeLJ>* pc) {
	for (auto & m : _vecOfMolecules) {
		pc->addParticle(m);
	}
}


void SingleCellIteratorTest::TearDown() {
}

TEST_F(SingleCellIteratorTest, testFullParticleCell) {
	FullParticleCell<MoleculeLJ> fpc;

	fillWithParticles(&fpc);

    SingleCellIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&fpc);
    int i = 0;
	for (; iter.isValid(); ++iter, ++i) {
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
		}
		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
	}
}

TEST_F(SingleCellIteratorTest, testRMMParticleCell) {
	RMMParticleCell<MoleculeLJ> fpc;

	fillWithParticles(&fpc);

    SingleCellIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&fpc);
    int i = 0;
	for (; iter.isValid(); ++iter, ++i) {
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
		}

//		IDs are not set by the RMM Cell yet!
//		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
	}
}
