/**
 * @file VerletListsTest.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

class VerletListsTest : public AutoPasTestBase, public ::testing::WithParamInterface<double> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto cellSizeFactor = static_cast<ParamType>(info.param);
      return "CellSizeFactor_" + std::to_string((int)cellSizeFactor);
    }
  };
};
