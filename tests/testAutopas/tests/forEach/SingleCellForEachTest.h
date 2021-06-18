/**
 * @file SingleCellIteratorTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

class SingleCellForEachTest : public AutoPasTestBase {
 public:
  SingleCellForEachTest() = default;

  void SetUp() override;

  void TearDown() override;

  ~SingleCellForEachTest() override = default;

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
