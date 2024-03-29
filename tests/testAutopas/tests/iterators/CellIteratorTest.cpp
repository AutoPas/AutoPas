/**
 * @file SingleCellIteratorTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "CellIteratorTest.h"

#include <gtest/gtest.h>

using namespace autopas;

void CellIteratorTest::SetUp() {
  for (int i = 0; i < 4; ++i) {
    std::array<double, 3> arr{};
    for (auto &a : arr) {
      a = static_cast<double>(i);
    }
    Molecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), 0);
    _vecOfMolecules.push_back(m);
  }
}

void CellIteratorTest::TearDown() {}

TEST_F(CellIteratorTest, testFullParticleCell) {
  FMCell fpc;

  fillWithParticles(&fpc);

  auto iter = fpc.begin();
  int i = 0;
  for (; iter != fpc.end(); ++iter, ++i) {
    for (int d = 0; d < 3; ++d) {
      ASSERT_DOUBLE_EQ(iter->getR()[d], _vecOfMolecules[i].getR()[d]);
    }
    ASSERT_EQ(iter->getID(), _vecOfMolecules[i].getID());
  }
}
