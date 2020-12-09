/**
 * @file SingleCellIteratorTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

namespace SingleCellIteratorTest {

class SingleCellIteratorTest : public AutoPasTestBase {
 public:
  SingleCellIteratorTest() = default;

  void SetUp() override;

  void TearDown() override;

  ~SingleCellIteratorTest() override = default;

  template <class Cell>
  void fillWithParticles(Cell *pc) {
    for (auto &m : _vecOfMolecules) {
      pc->addParticle(m);
    }
  }

 protected:
  // needs to be protected, because the test fixtures generate a derived class
  // for each unit test.
  std::vector<Molecule> _vecOfMolecules;
};

} // end namespace SingleCellIteratorTest
