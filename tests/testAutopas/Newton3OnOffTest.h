/**
 * @file Newton3OnOffTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPas.h"
#include "AutoPasTestBase.h"
#include "mocks/MockFunctor.h"
#include "sph/autopassph.h"
/**
 * Test to check if newton3 and non-newton3 work as expected
 */
class Newton3OnOffTest : public AutoPasTestBase {
 public:
  Newton3OnOffTest() : mockFunctor(), autoPas() {}

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }
  double getVerletSkin() const { return 0.0; }
  unsigned int getVerletRebuildFrequency() const { return 1; }

 protected:
  double fRand(double fMin, double fMax) const;

  std::array<double, 3> randomPosition(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) const;

  void fillContainerWithMolecules(unsigned long numMolecules,
                                  AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> &cont) const;

 public:
  MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> mockFunctor;
  AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autoPas;
};
