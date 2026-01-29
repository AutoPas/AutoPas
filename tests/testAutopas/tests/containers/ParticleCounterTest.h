/**
 * @file ParticleCounterTest.h
 * @author D. Martin
 * @date 10.09.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "testingHelpers/commonTypedefs.h"

class ParticleCounterTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::ContainerOption> {
 public:
  ParticleCounterTest() = default;

  ~ParticleCounterTest() override = default;

  static auto getTestParams();

 protected:
  static void calculateForces(autopas::ContainerOption containerOption, double cellSizeFactor);

  static constexpr std::array<double, 3> _boxMin{0, 0, 0};
  static constexpr std::array<double, 3> _boxMax{5, 5, 5};
  static constexpr double _cutoff{1.};
  static constexpr double _cellSizeFactor{1.};

  static constexpr unsigned int _numOwnedParticles = 3;
  static constexpr unsigned int _numHaloParticles = 4;
  static constexpr unsigned int _numDummyParticles = 1;

  static constexpr bool _useMortonIndex = true;
  static constexpr bool _preloadLJMixingPtr = true;
  static constexpr bool _useLiveId = true;
  static constexpr bool _reserveVLSizes = true;
};