/**
 * @file AutoPasVerletLikeInterfaceTest.h
 * @author seckler
 * @date 13.05.19
 */

#pragma once

#include <gtest/gtest.h>
#include "autopas/AutoPas.h"

class AutoPasInterfaceTest : public testing::Test,
                             public ::testing::WithParamInterface<autopas::ContainerOption> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string containerSelector_str(autopas::utils::StringUtils::to_string(inputTuple));
      return containerSelector_str;
    }
  };
};
