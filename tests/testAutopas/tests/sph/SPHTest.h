/**
 * @file SPHTest.h
 * @author seckler
 * @date 22.01.18
 */

#pragma once

#include <gtest/gtest.h>

#include <tuple>

#include "AutoPasTestBase.h"
#include "autopas/sph/autopassph.h"
#include "autopasTools/generators/RandomGenerator.h"

namespace SPHTest {

enum SPHFunctorType { density, hydro };

class SPHTest : public AutoPasTestBase,
                public ::testing::WithParamInterface<std::tuple<autopas::DataLayoutOption, SPHFunctorType>> {
  void SetUp() override{};

  void TearDown() override{};

 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [dataLayoutOption, sphFunctorType] = static_cast<ParamType>(info.param);
      std::string dataLayoutOptionStr(dataLayoutOption.to_string());
      std::string sphFunctorTypeStr(sphFunctorType == hydro ? "hydro" : "density");
      // replace all '-' with '_', otherwise the test name is invalid
      // std::replace(traversal.begin(), traversal.end(), '-', '_');
      std::string name = sphFunctorTypeStr + "Functor_" + dataLayoutOptionStr;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    }
  };
};

} // end namespace SPHTest
