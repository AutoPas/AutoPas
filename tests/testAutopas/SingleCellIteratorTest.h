/*
 * SingleCellIteratorTest.h
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#ifndef TESTS_TESTAUTOPAS_SINGLECELLITERATORTEST_H_
#define TESTS_TESTAUTOPAS_SINGLECELLITERATORTEST_H_

#include "autopas.h"
#include "gtest/gtest.h"

class SingleCellIteratorTest : public ::testing::Test {
public:
	SingleCellIteratorTest() {}

	void SetUp() override;

	void TearDown() override;

	~SingleCellIteratorTest() {}

	void testFullParticleCell();

	void testRMMParticleCell();

	void fillWithParticles(autopas::ParticleCell<autopas::MoleculeLJ> * pc);

protected:
	// needs to be protected, because the test fixtures generate a derived class for each unit test.
   std::vector<autopas::MoleculeLJ> _vecOfMolecules;
};


#endif /* TESTS_TESTAUTOPAS_SINGLECELLITERATORTEST_H_ */
