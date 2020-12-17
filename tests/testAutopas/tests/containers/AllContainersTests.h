/**
 * @file TestsAllContainers.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

using ParamType = autopas::ContainerOption;

class AllContainersTests : public AutoPasTestBase, public ::testing::WithParamInterface<ParamType> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      return info.param.to_string();
    };
    return paramToString;
  }

 protected:
  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {10, 10, 10};

  template <class ParticleType = autopas::Particle>
  auto getInitializedContainer() {
    auto containerOptionToTest = GetParam();
    double cutoff = 1;
    double skin = 0.2;
    double cellSizeFactor = 1;

    autopas::ContainerSelector<Particle> selector{boxMin, boxMax, cutoff};
    autopas::ContainerSelectorInfo selectorInfo{cellSizeFactor, skin, 32, autopas::LoadEstimatorOption::none};
    selector.selectContainer(containerOptionToTest, selectorInfo);
    return selector.getCurrentContainer();
  }

  /**
   * Adds three particles (lower corner, mid, and high corner) to the container
   * and also returns the particles for checks.
   * @tparam ParticleType
   * @param container
   * @return Vector of added particles
   */
  template <class ParticleType = autopas::Particle>
  std::vector<Particle> addParticlesMinMidMax(
      const std::shared_ptr<autopas::ParticleContainerInterface<Particle>> &container) {
    constexpr size_t numParticles = 3;
    std::vector<Particle> addedParticles;
    addedParticles.reserve(numParticles);

    // helper for Positions
    auto nearBoxMin = autopas::utils::ArrayMath::add(container->getBoxMin(), {0.1, 0.1, 0.1});
    auto nearBoxMax = autopas::utils::ArrayMath::sub(container->getBoxMax(), {0.1, 0.1, 0.1});
    auto nearBoxLength = autopas::utils::ArrayMath::sub(nearBoxMax, nearBoxMin);

    for (size_t i = 0; i < numParticles; ++i) {
      // place 3 particles:
      //   - lower corner
      //   - middle of domain
      //   - high corner
      // boxMin + i/2 * boxMax
      auto pos =
          autopas::utils::ArrayMath::add(nearBoxMin, autopas::utils::ArrayMath::mulScalar(nearBoxLength, i * 0.5));
      Particle particle(pos, {0., 0., 0.}, 1);
      container->addParticle(particle);
      addedParticles.push_back(particle);
    }
    return addedParticles;
  }

  void testUpdateContainerDeletesDummy(bool previouslyOwned);
};
