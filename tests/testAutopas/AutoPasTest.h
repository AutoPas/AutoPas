/**
 * @file AutoPasTest.h
 * @author seckler
 * @date 29.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPas.h"

class AutoPasTest : public testing::Test {
 public:
  AutoPasTest() {
    autoPas.init({0., 0., 0.}, {10., 10., 10.}, 1., 0, 1, {autopas::ContainerOptions::linkedCells},
                 {autopas::TraversalOptions::c08});
  }

  AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autoPas;
};
