/**
 * @file AoSvsSoATest.h
 * @author F.Gratl
 * @date 8.02.18
 */

#pragma once

#include <gtest/gtest.h>
#include <chrono>
#include "../../../../examples/md-flexible/ParticleClassLibrary.h"
#include "AutoPasTestBase.h"
#include "autopas/autopasIncludes.h"

class AoSvsSoATest : public AutoPasTestBase {
 public:
  void generateParticles(std::vector<autopas::Particle> *particles);
};
