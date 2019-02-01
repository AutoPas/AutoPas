/**
 * @file AutoPasTest.h
 * @author seckler
 * @date 29.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include "autopas/AutoPas.h"

class AutoPasTest : public testing::Test {
 public:
  AutoPasTest() {
    autoPas.init();
  }

  autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autoPas;
};
