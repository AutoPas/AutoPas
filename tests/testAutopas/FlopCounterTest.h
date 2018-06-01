/**
 * FlopCounterTest.h
 *
 *  Created on: 6/1/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include "AutoPasTestBase.h"
#include <AutoPas.h>

typedef autopas::Particle Particle;
typedef autopas::FullParticleCell<autopas::Particle> FPCell;

class FlopCounterTest : public AutoPasTestBase {
 public:
  FlopCounterTest() = default;

  ~FlopCounterTest() override = default;

  void test(autopas::DataLayoutOption dataLayoutOption);
};



