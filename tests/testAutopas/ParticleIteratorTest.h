/*
 * ParticleIteratorTest.h
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#ifndef TESTS_TESTAUTOPAS_PARTICLEITERATORTEST_H_
#define TESTS_TESTAUTOPAS_PARTICLEITERATORTEST_H_

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopasIncludes.h"

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

  std::vector<autopas::MoleculeLJ> _vecOfMolecules;
  unsigned long _currentIndex;
};

#endif /* TESTS_TESTAUTOPAS_PARTICLEITERATORTEST_H_ */
