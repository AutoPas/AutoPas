/**
 * @file RemainderTraversalTest3B.h
 * @author muehlhaeusser
 * @date 23.09.2023
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/GenerateValidConfigurations.h"

enum ParticleStorage {
  container,
  containerHalo,
  buffer,
  bufferHalo,
};

class RemainderTraversalTest3B
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<
          std::tuple<ParticleStorage, ParticleStorage, ParticleStorage, ContainerConfiguration,
                     autopas::DataLayoutOption, autopas::Newton3Option>> {
 public:
  RemainderTraversalTest3B() : numBuffers(autopas::autopas_get_max_threads()){};
  ~RemainderTraversalTest3B() override = default;

  size_t numBuffers;

  struct threeParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[choiceA, choiceB, choiceC, containerConfig, dataLayout, newton3] =
          static_cast<ParamType>(info.param);
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
      return enumToString(choiceA) + "_" + enumToString(choiceB) + "_" + enumToString(choiceC) + "_" +
             containerConfig.toShortString() + "_" + dataLayout.to_string() + "_N3_" + newton3.to_string();
    }
  };
};

/**
 * Fixture for testing the individual steps of the triwise remainder traversal directly, where each test places exactly
 * three particles in specific storage locations (container, container halo, particle buffers, halo particle buffers).
 *
 * Tests are parameterized over all combinations of container configuration, data layout, and Newton3 option for which a
 * valid triwise configuration exists.
 */
class RemainderTraversalDirectTest3B
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<
          std::tuple<ContainerConfiguration, autopas::DataLayoutOption, autopas::Newton3Option>> {
 public:
  RemainderTraversalDirectTest3B() : numBuffers(autopas::autopas_get_max_threads()){};
  ~RemainderTraversalDirectTest3B() override = default;

  /**
   * Number of particle buffers, i.e. one per OpenMP thread.
   */
  size_t numBuffers;

  /**
   * Generates a test name of the form "<container>_csf_<csf>_<dataLayout>_N3_<newton3>" from the test parameters.
   */
  struct paramToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[containerConfig, dataLayout, newton3] = static_cast<ParamType>(info.param);
      return containerConfig.toShortString() + "_" + dataLayout.to_string() + "_N3_" + newton3.to_string();
    }
  };
};
