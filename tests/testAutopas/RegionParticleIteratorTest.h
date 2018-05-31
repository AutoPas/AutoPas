/**
 * @file RegionParticleIteratorTest.h
 * @author seckler
 * @date 03.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopasIncludes.h"

class TouchableParticle : public autopas::Particle {
 public:
  TouchableParticle(std::array<double, 3> pos, unsigned long id)
      : autopas::Particle(pos, {0, 0, 0}, id), _numTouched(0){};
  void touch() { _numTouched++; }
  unsigned int getNumTouched() { return _numTouched; }

 private:
  unsigned int _numTouched;
};

class RegionParticleIteratorTest : public AutoPasTestBase {
 public:
  RegionParticleIteratorTest()
      : _boxMin{0., 0., 0.},
        _boxMax{5., 5., 5.},
        _regionMin{1., 1., 1.},
        _regionMax{3., 3., 3.},
        _cutoff{.9} {}

  void SetUp() override{};

  void TearDown() override{};

  ~RegionParticleIteratorTest() override = default;

 protected:
  // needs to be protected, because the test fixtures generate a derived class
  // for each unit test.

  std::array<double, 3> _boxMin, _boxMax;

  std::array<double, 3> _regionMin, _regionMax;

  double _cutoff;
};
