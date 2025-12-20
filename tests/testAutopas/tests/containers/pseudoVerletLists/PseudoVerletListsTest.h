/**
* @file PseudoVerletListsTest.h
 * @author Lars Doll
 * @date 20.12.2025
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/containers/pseudoVerletLists/PseudoVerletLists.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"

class PseudoVerletListsTest : public AutoPasTestBase, public ::testing::WithParamInterface<double> {
public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto cellSizeFactor = static_cast<ParamType>(info.param);
      return "CellSizeFactor_" + std::to_string((int)cellSizeFactor);
    }
  };
};
