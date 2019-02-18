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
    autoPas.setBoxMin({0., 0., 0.});
    autoPas.setBoxMax({10., 10., 10.});
    autoPas.setCutoff(1.);
    autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
    autoPas.setAllowedTraversals({autopas::TraversalOption::c08});
    autoPas.init();
  }

  autopas::AutoPas<autopas::Particle, autopas::FullParticleCell<autopas::Particle>> autoPas;
};
