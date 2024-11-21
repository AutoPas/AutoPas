

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayMath.h"
#include "src/TypeDefinitions.h"
#include "src/zonalMethods/FullShell.h"
#include "src/zonalMethods/region/RectRegion.h"

using namespace autopas::utils::ArrayMath::literals;

/**
 * Test class for the Zone classes
 */
class FullShellTest : public AutoPasTestBase, public FullShell {
 public:
  using AutoPasType = autopas::AutoPas<ParticleType>;

  /**
   * Constructor.
   */
  FullShellTest()
      : AutoPasTestBase(),
        FullShell(_cutoff, _verletSkinWidth, 0, _homeBoxRegion, _homeBoxRegion)

  {}

 protected:
  // initialize autopas container
  void initContainer(AutoPasType &autopas, std::vector<ParticleType> particles);

  AutoPasType _autopas;

  constexpr static double _cutoff = 0.9;
  constexpr static double _verletSkinWidth = 0.1;
  constexpr static std::array<double, 3> _boxMin = {0, 0, 0};
  constexpr static std::array<double, 3> _boxMax = {10, 10, 10};
  inline const static RectRegion _homeBoxRegion{_boxMin, _boxMax - _boxMin};
};
