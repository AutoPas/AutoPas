/**
 * @file RemainderTraversalTest.h
 * @author F. Gratl
 * @date 28.11.2022
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

class RemainderTraversalTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<std::tuple<ParticleStorage, ParticleStorage, ContainerConfiguration,
                                                      autopas::DataLayoutOption, autopas::Newton3Option>> {
 public:
  RemainderTraversalTest() : numBuffers(autopas::autopas_get_max_threads()){};
  ~RemainderTraversalTest() override = default;

  size_t numBuffers;

  struct twoParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &[choiceA, choiceB, containerConfig, dataLayout, newton3] = static_cast<ParamType>(info.param);
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
      return enumToString(choiceA) + "_" + enumToString(choiceB) + "_" + containerConfig.toShortString() + "_" +
             dataLayout.to_string() + "_N3_" + newton3.to_string();
    }
  };
};

/**
 * Fixture for testing the individual steps of the remainder traversal directly, where each test places exactly two
 * particles in specific storage locations (container, container halo, particle buffers, halo particle buffers).
 *
 * Tests are parameterized over all combinations of container configuration, data layout, and Newton3 option for which a
 * valid configuration exists.
 */
class RemainderTraversalDirectTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<
          std::tuple<ContainerConfiguration, autopas::DataLayoutOption, autopas::Newton3Option>> {
 public:
  RemainderTraversalDirectTest() : numBuffers(autopas::autopas_get_max_threads()){};
  ~RemainderTraversalDirectTest() override = default;

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
