/**
 * @file Newton3OnOffTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/sph/autopassph.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/RandomGenerator.h"

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
 public:
  MockFunctor<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> mockFunctor;
  autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autoPas;
};
