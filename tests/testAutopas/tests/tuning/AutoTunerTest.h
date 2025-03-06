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
      autopas::ContainerOption::directSum,      _cellSizeFactor,
      autopas::TraversalOption::ds_sequential,  autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::enabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confDs_seq_noN3{
      autopas::ContainerOption::directSum,      _cellSizeFactor,
      autopas::TraversalOption::ds_sequential,  autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confLc_c01_noN3{
      autopas::ContainerOption::linkedCells,    _cellSizeFactor,
      autopas::TraversalOption::lc_c01,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confLc_c18_noN3{
      autopas::ContainerOption::linkedCells,    _cellSizeFactor,
      autopas::TraversalOption::lc_c18,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confLc_c18_N3{
      autopas::ContainerOption::linkedCells,    _cellSizeFactor,
      autopas::TraversalOption::lc_c18,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::enabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confLc_c08_N3{
      autopas::ContainerOption::linkedCells,    _cellSizeFactor,
      autopas::TraversalOption::lc_c08,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::enabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confLc_c08_noN3{
      autopas::ContainerOption::linkedCells,    _cellSizeFactor,
      autopas::TraversalOption::lc_c08,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  // Triwise configs:
  const autopas::Configuration _confLc_c01_3b_noN3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c01,        autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise, autopas::VectorizationPatternOption::p1xVec};
  const autopas::Configuration _confDs_3b_N3{
      autopas::ContainerOption::directSum,     _cellSizeFactor,
      autopas::TraversalOption::ds_sequential, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,          autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::triwise, autopas::VectorizationPatternOption::p1xVec};
};
