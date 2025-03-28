/**
 * @file AutoPasInterfaceTest.h
 * @author seckler
 * @date 13.05.19
 */

#pragma once

#include <gtest/gtest.h>

#include <tuple>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/utils/ArrayMath.h"

class AutoPasInterfaceTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::Configuration> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto conf = static_cast<ParamType>(info.param);
      std::string str;
      str += conf.container.to_string() + "_";
      str += conf.traversal.to_string() + "_";
      str += conf.loadEstimator.to_string() + "_";
      str += conf.dataLayout.to_string() + "_";
      str += "N3" + conf.newton3.to_string() + "_";
      str += std::string{"cellSizeFactor"} + std::to_string(conf.cellSizeFactor);
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };
};

class AutoPasInterface1ContainersTest : public testing::Test,
                                        public ::testing::WithParamInterface<autopas::ContainerOption> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      const auto &option = static_cast<ParamType>(info.param);
      return option.to_string();
    }
  };
};

class AutoPasInterface2ContainersTest
    : public testing::Test,
      public ::testing::WithParamInterface<std::tuple<autopas::ContainerOption, autopas::ContainerOption>> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);
      std::string fromStr(std::get<0>(inputTuple).to_string());
      std::string toStr(std::get<1>(inputTuple).to_string());
      // replace all '-' with '_', otherwise the test name is invalid
      // std::replace(traversal.begin(), traversal.end(), '-', '_');
      return fromStr + "And" + toStr;
    }
  };
};

static inline auto getTestableContainerOptions() { return autopas::ContainerOption::getAllOptions(); }

/**
 * Adds three particles (lower corner, mid, and high corner) to the container
 * and also returns copies of the particles for checks.
 * @tparam ParticleType
 * @param container
 * @return Vector of added particles
 */
template <class Container>
std::vector<typename Container::ParticleType> addParticlesMinMidMax(Container &container) {
  using namespace autopas::utils::ArrayMath::literals;

  constexpr size_t numParticles = 3;
  std::vector<typename Container::ParticleType> addedParticles;
  addedParticles.reserve(numParticles);

  // helper for Positions
  auto nearBoxMin = autopas::utils::ArrayMath::add(container.getBoxMin(), {0.1, 0.1, 0.1});
  auto nearBoxMax = autopas::utils::ArrayMath::sub(container.getBoxMax(), {0.1, 0.1, 0.1});
  auto nearBoxLength = nearBoxMax - nearBoxMin;

  for (size_t i = 0; i < numParticles; ++i) {
    // place 3 particles:
    //   - lower corner
    //   - middle of domain
    //   - high corner

    // factor for relative positioning in the domain of the i-th particle
    constexpr double scalingFactor = 1. / (numParticles - 1);
    // boxMin + i/scalingFactor * boxMax
    auto pos = nearBoxMin + (nearBoxLength * (i * scalingFactor));
    typename Container::ParticleType particle(pos, {0., 0., 0.}, 1);
    container.addParticle(particle);
    addedParticles.push_back(particle);
  }
  return addedParticles;
}