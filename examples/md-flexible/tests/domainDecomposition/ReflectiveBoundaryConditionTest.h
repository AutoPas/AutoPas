/**
 * @file ReflectiveBoundaryConditionTest.h
 * @author F. Gratl
 * @date 26.04.22
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/Quaternion.h"
#include "autopas/utils/WrapMPI.h"
/**
 * Parameterized test case for reflective boundary conditions in RegularGridDecomposition
 * Quaternion parameter is not used for single-site compilations.
 */
class ReflectiveBoundaryConditionTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple</*position*/ std::array<double, 3>,
                                                      /*quaternion*/ std::array<double, 4>>> {
 public:
  /**
   * Constructor.
   */
  ReflectiveBoundaryConditionTest() = default;

  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &inputTuple = static_cast<ParamType>(info.param);
      const auto &[pos, vel, isReflected] = inputTuple;
      std::stringstream ss;
      using autopas::utils::ArrayUtils::to_string;
      ss << "pos" << to_string(pos, "_", {"", ""}) << "_vel" << to_string(vel, "_", {"", ""}) << "_dimReflects"
         << to_string(isReflected, "_", {"", ""});
      auto res = ss.str();
      std::replace(res.begin(), res.end(), '-', '_');
      std::replace(res.begin(), res.end(), '.', '_');
      return res;
    }
  };
};