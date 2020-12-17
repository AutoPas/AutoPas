/**
 * @file ContainerSelectorTest.h
 * @author F. Gratl
 * @date 22.06.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

class ContainerSelectorTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::ContainerOption> {
 public:
  ContainerSelectorTest() = default;
  ~ContainerSelectorTest() override = default;

  struct oneParamToString {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      // tuple of ContainerOption
      const auto &option = static_cast<ParamType>(info.param);
      return option.to_string();
    }
  };

 protected:
  const std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0.1;

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
};
