/**
 * @file AutoTunerTest.h
 * @author F. Gratl
 * @date 8/10/18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class AutoTunerTest : public AutoPasTestBase {
 public:
  AutoTunerTest() = default;
  ~AutoTunerTest() override = default;

  const double _cellSizeFactor{1.};

  // configurations used throughout various tests
  const autopas::Configuration _confDs_seq_N3{
      autopas::ContainerOption::directSum,     _cellSizeFactor,
      autopas::TraversalOption::ds_sequential, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::enabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confDs_seq_noN3{
      autopas::ContainerOption::directSum,     _cellSizeFactor,
      autopas::TraversalOption::ds_sequential, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confLc_c01_noN3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c01,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confLc_c18_noN3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c18,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confLc_c18_N3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c18,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::enabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confLc_c08_N3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c08,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::enabled,
      autopas::InteractionTypeOption::pairwise};
  const autopas::Configuration _confLc_c08_noN3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c08,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise};
  // Triwise configs:
  const autopas::Configuration _confLc_c01_3b_noN3{
      autopas::ContainerOption::linkedCells,  _cellSizeFactor,
      autopas::TraversalOption::lc_c01,       autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,         autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise};
  const autopas::Configuration _confDs_3b_N3{
      autopas::ContainerOption::directSum,     _cellSizeFactor,
      autopas::TraversalOption::ds_sequential, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise};

 protected:
  template <typename Functor>
  void testAllConfigsOfType(Functor &functor, size_t numExpectedConfigs,
                            const autopas::InteractionTypeOption &interactionType);

  /**
   * Tests that ending a tuning phase with a reject configuration is handled correctly.
   *
   * Mimics an AutoTuner trialling four configurations. The second and the final configurations will be rejected. After
   * the tuning phase is completed, the test should get the correct configuration upon calling `getNextConfig`
   *
   * @param rejectIndefinitely Run the test, parsing this rejectIndefinitely value to rejectConfig
   */
  void testEndingTuningPhaseWithRejectedConfig(bool rejectIndefinitely) const;
};
