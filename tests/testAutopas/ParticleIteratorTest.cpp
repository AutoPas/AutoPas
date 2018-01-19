/*
 * ParticleIteratorTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#include "ParticleIteratorTest.h"
#include "autopas.h"
#include "gtest/gtest.h"

using namespace autopas;

void ParticleIteratorTest::SetUp() {
	for (int i = 0; i < 20; ++i) {
		std::array<double, 3> arr;
		for (auto & a : arr) {
			a = static_cast<double>(i);
		}
		MoleculeLJ m(arr, { 0., 0., 0. }, static_cast<unsigned long>(i));
		_vecOfMolecules.push_back(m);
	}

	_currentIndex = 0;
}

void ParticleIteratorTest::TearDown() {

}

void ParticleIteratorTest::fillWithParticles(
		autopas::ParticleCell<autopas::MoleculeLJ>* pc) {

	// insert four particles
	for (int i = _currentIndex; i < _currentIndex + 4; ++i) {
		pc->addParticle(_vecOfMolecules.at(i));
	}
	_currentIndex += 4;
}

TEST_F(ParticleIteratorTest, testFullIterator_EFEFFEEFEF) {
	// Empty Full Empty Full Full Empty Empty Full Empty Full
	std::vector<FullParticleCell<MoleculeLJ>> data(10);

	for (auto i : {1, 3, 4, 7, 9}) {
		fillWithParticles(&data.at(i));
	}

    ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
	for (; iter.isValid(); ++iter, ++i) {
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
		}
		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
	}

}

TEST_F(ParticleIteratorTest, testFullIterator_FEFEEFFEFE) {
	// Full Empty Full Empty Empty Full Full Empty Full Empty
	std::vector<FullParticleCell<MoleculeLJ>> data(10);

	for (auto i : {0, 2, 5, 6, 8}) {
		fillWithParticles(&data.at(i));
	}

    ParticleIterator<MoleculeLJ, FullParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
	for (; iter.isValid(); ++iter, ++i) {
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
		}
		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
	}
}

TEST_F(ParticleIteratorTest, testRMMIterator_EFEFFEEFEF) {
	// Empty Full Empty Full Full Empty Empty Full Empty Full
	std::vector<RMMParticleCell<MoleculeLJ>> data(10);

	for (auto i : {1, 3, 4, 7, 9}) {
		fillWithParticles(&data.at(i));
	}

    ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
	for (; iter.isValid(); ++iter, ++i) {
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
		}
//		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
	}

}

TEST_F(ParticleIteratorTest, testRMMIterator_FEFEEFFEFE) {
	// Full Empty Full Empty Empty Full Full Empty Full Empty
	std::vector<RMMParticleCell<MoleculeLJ>> data(10);

	for (auto i : {0, 2, 5, 6, 8}) {
		fillWithParticles(&data.at(i));
	}

    ParticleIterator<MoleculeLJ, RMMParticleCell<MoleculeLJ>> iter(&data);
    int i = 0;
	for (; iter.isValid(); ++iter, ++i) {
		for (int d = 0; d < 3; ++d) {
			ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
		}
//		ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
	}
}
