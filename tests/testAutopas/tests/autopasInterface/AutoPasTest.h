/**
 * @file AutoPasTest.h
 * @author seckler
 * @date 29.05.18
 */

#pragma once

#include <gtest/gtest.h>

#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

namespace AutoPasTest {

class AutoPasTest : public testing::Test {
 public:
  AutoPasTest() {
    autoPas.setBoxMin({0., 0., 0.});
    autoPas.setBoxMax({10., 10., 10.});
    autoPas.setCutoff(1.);
    autoPas.init();
  }

 protected:
  void expectedParticles(size_t expectedOwned, size_t expectedHalo);

  autopas::AutoPas<Particle> autoPas;
};

} // end namespace AutoPasTest
