/**
 * @file AoSvsSoATest.h
 * @author F.Gratl
 * @date 8.02.18
 */

#pragma once

#include <gtest/gtest.h>
#include <chrono>
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"
#include "autopas/molecularDynamics/ParticleClassLibrary.h"

class AoSvsSoATest : public AutoPasTestBase {
 public:
  void generateParticles(std::vector<autopas::Particle> *particles);
};
