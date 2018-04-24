/**
 * @file RegionParticleIteratorTest.h
 * @author seckler
 * @date 03.04.18
 */

#pragma once

#include "autopasIncludes.h"
#include "gtest/gtest.h"
#include "AutoPasTest.h"

class TouchableParticle : public autopas::Particle {
 public:
  TouchableParticle(std::array<double, 3> pos, unsigned long id)
      : autopas::Particle(pos, {0, 0, 0}, id), _numTouched(0){};
  void touch() { _numTouched++; }
  unsigned int getNumTouched() { return _numTouched; }

 private:
  unsigned int _numTouched;
};

class RegionParticleIteratorTest : public AutoPasTest {
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

  template <class Container>
  void fillWithParticles(Container &container) {
    srand(42);  // fixed seedpoint
    int numParticles = 100;

    for (int i = 0; i < numParticles; ++i) {
      auto id = static_cast<unsigned long>(i);
      TouchableParticle particle(
          randomPosition(container.getBoxMin(), container.getBoxMax()), id);
      container.addParticle(particle);
    }
  }

  double fRand(double fMin, double fMax) {
    double f = static_cast<double>(rand()) / RAND_MAX;
    return fMin + f * (fMax - fMin);
  }

  std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin,
                                       const std::array<double, 3> &boxMax) {
    std::array<double, 3> r{0, 0, 0};
    for (int d = 0; d < 3; ++d) {
      r[d] = fRand(boxMin[d], boxMax[d]);
    }
    return r;
  }

  std::array<double, 3> _boxMin, _boxMax;

  std::array<double, 3> _regionMin, _regionMax;

  double _cutoff;
};
