/**
 * @file TimeDiscretizationTest.h
 * @author N. Fottner
 * @date 05/22/19.
 */

#pragma once

#include <gtest/gtest.h>

#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

class TimeDiscretizationTest : public AutoPasTestBase {
 public:
  TimeDiscretizationTest() : AutoPasTestBase(), _particlePropertiesLibrary(1) {
    _particlePropertiesLibrary.addType(0, 1, 1, 1);
    _particlePropertiesLibrary.calculateMixingCoefficients();
  }
  void fillWithParticlesAndInit(autopas::AutoPas<Molecule, FMCell> &autopas);

 protected:
  ParticlePropertiesLibrary<double, size_t> _particlePropertiesLibrary;
};