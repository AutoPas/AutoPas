/*
 * CellBlock3DTest.cpp
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#include "CellBlock3DTest.h"

void CellBlock3DTest::SetUp() {
	// init
	std::array<double, 3> bMin({0.0, 0.0, 0.0});
	std::array<double, 3> bMax({10.0, 10.0, 10.0});
	double intLength;

	intLength = 10.0;
	_cells_1x1x1.init(&_vec1, bMin, bMax, intLength);

	intLength = 5.0;
	_cells_2x2x2.init(&_vec2, bMin, bMax, intLength);

	intLength = 3.0;
	_cells_3x3x3.init(&_vec3, bMin, bMax, intLength);
}

TEST_F(CellBlock3DTest, test1x1x1) {
	std::array<double, 3> start = { -5., -5., -5. }, dr = { 10.0, 10.0, 10.0 };
	std::array<int, 3> numParts = { 3, 3, 3 };
	auto mesh = getMesh(start, dr, numParts);

	int counter = 0;
	for (auto & m : mesh) {
		int index = _cells_1x1x1.get1DIndexOfPosition(m);
		ASSERT_EQ(index, counter);
		++counter;
	}
}

TEST_F(CellBlock3DTest, test2x2x2) {
	std::array<double, 3> start = { -2.5, -2.5, -2.5 }, dr = { 5.0, 5.0, 5.0 };
	std::array<int, 3> numParts = { 4, 4, 4 };
	auto mesh = getMesh(start, dr, numParts);

	int counter = 0;
	for (auto & m : mesh) {
		int index = _cells_2x2x2.get1DIndexOfPosition(m);
		ASSERT_EQ(index, counter);
		++counter;
	}
}

TEST_F(CellBlock3DTest, test3x3x3) {
	std::array<double, 3> start = { -1.6, -1.6, -1.6 }, dr = { 3.3, 3.3, 3.3 };
	std::array<int, 3> numParts = { 5, 5, 5 };
	auto mesh = getMesh(start, dr, numParts);

	int counter = 0;
	for (auto & m : mesh) {
		int index = _cells_3x3x3.get1DIndexOfPosition(m);
		ASSERT_EQ(index, counter);
		++counter;
	}
}

std::vector<std::array<double, 3> > CellBlock3DTest::getMesh(
		std::array<double, 3> start, std::array<double, 3> dr,
		std::array<int, 3> numParts) const {
	std::vector<std::array<double, 3> > ret;

	for (int z = 0; z < numParts[2]; ++z) {
		for (int y = 0; y < numParts[1]; ++y) {
			for (int x = 0; x < numParts[0]; ++x) {
				std::array<double, 3> pos;
				pos[0] = start[0] + x * dr[0];
				pos[1] = start[1] + y * dr[1];
				pos[2] = start[2] + z * dr[2];
				ret.push_back(pos);
			}
		}
	}
	return ret;
}
