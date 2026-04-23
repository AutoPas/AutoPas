/**
 * @file TuningManagerTest.h
 * @author muehlhaeusser
 * @date 18.04.2026
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/tuning/Configuration.h"

class TuningManagerTest : public AutoPasTestBase {
 public:
  TuningManagerTest() = default;
  ~TuningManagerTest() override = default;

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
  const autopas::Configuration _confLc_c01_noN3{
      autopas::ContainerOption::linkedCells,   _cellSizeFactor,
      autopas::TraversalOption::lc_c01,        autopas::LoadEstimatorOption::none,
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
};
