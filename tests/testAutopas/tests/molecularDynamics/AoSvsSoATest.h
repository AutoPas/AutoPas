/**
 * @file AoSvsSoATest.h
 * @author F.Gratl
 * @date 8.02.18
 */

#pragma once

#include <gtest/gtest.h>

#include <chrono>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

namespace AoSvsSoATest {

class AoSvsSoATest : public AutoPasTestBase {
 public:
  void generateParticles(std::vector<Molecule> *particles);
};

} // end namespace AoSvsSoATest
