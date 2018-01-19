/*
 * CellBlock3DTest.h
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#ifndef TESTS_TESTAUTOPAS_CELLBLOCK3DTEST_H_
#define TESTS_TESTAUTOPAS_CELLBLOCK3DTEST_H_

#include "autopas.h"
#include "gtest/gtest.h"

class CellBlock3DTest : public ::testing::Test {
public:
	CellBlock3DTest() {}

	void SetUp() override;

	~CellBlock3DTest() {}

	void test1x1x1();

	void test2x2x2();

	void test3x3x3();

protected:
	std::vector<std::array<double, 3>> getMesh(std::array<double, 3> start,
			std::array<double, 3> dx, std::array<int, 3> numParts) const;

	std::vector<autopas::FullParticleCell<autopas::MoleculeLJ>> _vec1, _vec2, _vec3;
	autopas::CellBlock3D<autopas::FullParticleCell<autopas::MoleculeLJ>>
		_cells_1x1x1, _cells_2x2x2, _cells_3x3x3;
};


#endif /* TESTS_TESTAUTOPAS_CELLBLOCK3DTEST_H_ */
