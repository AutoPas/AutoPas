/**
 * @file ParticleViewTest.h
 * @author lgaertner
 * @date 01.12.2021
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/KokkosParticleCell.h"
#include "autopas/kokkosContainers/ParticleView.h"
#include "testingHelpers/commonTypedefs.h"

class ParticleViewTest : public AutoPasTestBase {
 public:
  ParticleViewTest();
};
