/**
 * @file PairwiseVerletListsTest.h
 * @author tirgendetwas
 * @date 22.11.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/Particle.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"
#include "tests/containers/verletListsCellBased/verletLists/VerletListsTest.h"

class PairwiseVerletListsTest : public AutoPasTestBase, public ::testing::WithParamInterface<std::pair<double, bool>> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto pairParams = static_cast<ParamType>(info.param);
      return "CellSizeFactor_" + std::to_string((int)pairParams.first) + "_useNewton3_" +
             boolToString(pairParams.second);
    }

    std::string boolToString(bool n3) const {
      if (n3 == true) {
        return "true";
      } else {
        return "false";
      }
    }
  };
};