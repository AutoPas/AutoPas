/**
 * @file RemainderTraversalTest3B.h
 * @author muehlhaeusser
 * @date 23.09.2023
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/WrapOpenMP.h"

enum ParticleStorage {
  container,
  containerHalo,
  buffer,
  bufferHalo,
};

class RemainderTraversalTest3B
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<ParticleStorage, ParticleStorage, ParticleStorage>> {
 public:
  RemainderTraversalTest3B() : numBuffers(autopas::autopas_get_max_threads()){};
  ~RemainderTraversalTest3B() override = default;

  size_t numBuffers;

  struct threeParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[choiceA, choiceB, choiceC] = static_cast<ParamType>(info.param);
      auto enumToString = [](const auto &e) -> std::string {
        switch (e) {
          case ParticleStorage::container:
            return "container";
          case ParticleStorage::containerHalo:
            return "containerHalo";
          case ParticleStorage::buffer:
            return "buffer";
          case ParticleStorage::bufferHalo:
            return "bufferHalo";
          default:
            return "unknown";
        }
      };
      return enumToString(choiceA) + "_" + enumToString(choiceB) + "_" + enumToString(choiceC);
    }
  };
};
