/**
 * @file FmmTreeTest.h
 * @author Joachim Marin
 * @date 5.11.2019
 */

#pragma once

#include <gtest/gtest.h>
#include <array>
#include <iostream>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/particles/Particle.h"
#include "testingHelpers/RandomGenerator.h"

using AutoPasCont = autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>>;

class FmmTreeTest : public AutoPasTestBase {
 public:
  FmmTreeTest() = default;
};
