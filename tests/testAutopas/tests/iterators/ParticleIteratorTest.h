/**
 * @file ParticleIteratorTest.h
 * @author tchipev
 * @date 19.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/commonTypedefs.h"

namespace ParticleIteratorTest {

class ParticleIteratorTest : public AutoPasTestBase {
 public:
  ParticleIteratorTest() : _currentIndex(0ul) {}

  void SetUp() override;

  void TearDown() override;

  ~ParticleIteratorTest() override = default;

  template <class CellType>
  void fillWithParticles(CellType *pc) {
    // insert four particles
    for (unsigned long i = _currentIndex; i < _currentIndex + 4; ++i) {
      pc->addParticle(_vecOfMolecules.at(i));
    }
    _currentIndex += 4;
  }

 protected:
  // needs to be protected, because the test fixtures generate a derived class
  // for each unit test.

  std::vector<Molecule> _vecOfMolecules;
  unsigned long _currentIndex;
};

} // end namespace ParticleIteratorTest
