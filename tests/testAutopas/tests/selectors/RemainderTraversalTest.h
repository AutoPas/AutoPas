/**
 * @file RemainderTraversalTest.h
 * @author F. Gratl
 * @date 28.11.2022
 */

#pragma once

#include "AutoPasTestBase.h"

enum ParticleStorage {
  container,
  containerHalo,
  buffer,
  bufferHalo,
};

class RemainderTraversalTest : public AutoPasTestBase,
                               public ::testing::WithParamInterface<std::tuple<ParticleStorage, ParticleStorage>> {
 public:
  RemainderTraversalTest() = default;
  ~RemainderTraversalTest() override = default;

  struct twoParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[choiceA, choiceB] = static_cast<ParamType>(info.param);
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
      return enumToString(choiceA) + "_" + enumToString(choiceB);
    }
  };
};
