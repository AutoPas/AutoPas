/**
 * @file PairwiseVerletListsTest.h
 * @author tirgendetwas
 * @date 22.11.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "mocks/MockPairwiseFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class PairwiseVerletListsTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<double, bool, autopas::VerletListsCellsHelpers::VLCBuildType>> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto params = static_cast<ParamType>(info.param);
      return "CellSizeFactor_" + std::to_string((int)std::get<0>(params)) + "_useNewton3_" +
             boolToString(std::get<1>(params)) + "_buildType_" + buildTypeToString(std::get<2>(params));
    }

    std::string boolToString(bool n3) const {
      if (n3 == true) {
        return "true";
      } else {
        return "false";
      }
    }

    std::string buildTypeToString(autopas::VerletListsCellsHelpers::VLCBuildType buildType) const {
      if (buildType == autopas::VerletListsCellsHelpers::VLCBuildType::soaBuild) {
        return "buildSoA";
      } else {
        return "buildAoS";
      }
    }
  };
};