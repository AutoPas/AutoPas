/*
 * ParticleIteratorTest.h
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#ifndef TESTS_TESTAUTOPAS_PARTICLEITERATORTEST_H_
#define TESTS_TESTAUTOPAS_PARTICLEITERATORTEST_H_

#include "autopas.h"
#include "gtest/gtest.h"

class ParticleIteratorTest: public testing::Test {
public:
	ParticleIteratorTest() {}

	void SetUp() override;

	void TearDown() override;

	~ParticleIteratorTest() {}

	void fillWithParticles(autopas::ParticleCell<autopas::MoleculeLJ> * pc);

protected:
	// needs to be protected, because the test fixtures generate a derived class for each unit test.

	std::vector<autopas::MoleculeLJ> _vecOfMolecules;
	int _currentIndex;
};

#endif /* TESTS_TESTAUTOPAS_PARTICLEITERATORTEST_H_ */
