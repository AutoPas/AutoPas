/**
 * @file CellIteratorTest.h
 * @author F. Gratl
 * @date 12.01.23
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * There used to be a dedicated CellIterator.
 * In https://github.com/AutoPas/AutoPas/pull/712 this was removed so this test now only
 * tests if FullParticleCell provides good interfaces to std::vector::iterator.
 */
class CellIteratorTest : public AutoPasTestBase {
 public:
  CellIteratorTest() = default;

  void SetUp() override;

  void TearDown() override;

  ~CellIteratorTest() override = default;

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
