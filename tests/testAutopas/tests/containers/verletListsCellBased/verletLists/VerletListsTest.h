/**
 * @file VerletListsTest.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "generators/src/UniformGenerator.h"
#include "mocks/MockPairwiseFunctor.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

class VerletListsTest : public AutoPasTestBase, public ::testing::WithParamInterface<std::tuple<double, bool>> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [cellSizeFactor, newton3] = static_cast<ParamType>(info.param);
      std::string newton3Str = newton3 ? "N3" : "noN3";
      return "CellSizeFactor_" + std::to_string((int)cellSizeFactor) + "_" + newton3Str;
    }
  };
};
