/**
 * @file RegionParticleIteratorTest.h
 * @author seckler
 * @date 03.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/TouchableParticle.h"

class RegionParticleIteratorTest : public AutoPasTestBase {
 public:
  using LCTouch = autopas::LinkedCells<TouchableParticle>;

  RegionParticleIteratorTest()
      : _boxMin{0., 0., 0.}, _boxMax{5., 5., 5.}, _regionMin{1., 1., 1.}, _regionMax{3., 3., 3.}, _cutoff{.9} {}

  void SetUp() override{};

  void TearDown() override{};

  ~RegionParticleIteratorTest() override = default;

  /**
   * Checks if all particles in the given region of a Lined Cells container are touched exactly once and provides debug
   * output.
   * @param lcContainer
   * @param regionMin
   * @param regionMax
   */
  static void checkTouches(LCTouch &lcContainer, std::array<double, 3> &regionMin, std::array<double, 3> &regionMax);

 protected:
  void testLinkedCellsRegionParticleIteratorBehaviorOwned();

  void testLinkedCellsRegionParticleIteratorBehaviorHalo();

  // needs to be protected, because the test fixtures generate a derived class
  // for each unit test.

  std::array<double, 3> _boxMin, _boxMax;

  std::array<double, 3> _regionMin, _regionMax;

  double _cutoff;
};
