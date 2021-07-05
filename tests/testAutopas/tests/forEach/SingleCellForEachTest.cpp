/**
 * @file SingleCellIteratorTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "SingleCellForEachTest.h"

#include <gtest/gtest.h>

using namespace autopas;

void SingleCellForEachTest::SetUp() {
  for (int i = 0; i < 4; ++i) {
    std::array<double, 3> arr{};
    for (auto &a : arr) {
      a = static_cast<double>(i);
    }
    Molecule m(arr, {0., 0., 0.}, static_cast<unsigned long>(i), 0);
    _vecOfMolecules.push_back(m);
  }
}

void SingleCellForEachTest::TearDown() {}

TEST_F(SingleCellForEachTest, testFullParticleCell) {
  FMCell fpc;

  fillWithParticles(&fpc);

  auto iter = fpc.begin();
  int i = 0;

  auto forEachLambda = [&](auto p) {
    for (int d = 0; d < 3; ++d) {
      ASSERT_DOUBLE_EQ(p.getR()[d], _vecOfMolecules[i].getR()[d]);
    }
    ASSERT_EQ(p.getID(), _vecOfMolecules[i].getID());
    ++i;  // TODO lgaertner: wtf NO, but how else? maybe _vecOfMolecules.containsID(p.getID())
  };

  fpc.forEach(forEachLambda);
}
