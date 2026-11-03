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

  // Pairwise Linked Cells Config
  const double _cellSizeFactor{1.};
  const autopas::Configuration _confLc_c01_noN3{
      autopas::ContainerOption::linkedCells,    _cellSizeFactor,
      autopas::TraversalOption::lc_c01,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::NA};

  // Pairwise Verlet Lists Config
  const autopas::Configuration _confVl_list_iteration{autopas::ContainerOption::verletLists,
                                                      _cellSizeFactor,
                                                      autopas::TraversalOption::vl_list_iteration,
                                                      autopas::LoadEstimatorOption::none,
                                                      autopas::DataLayoutOption::aos,
                                                      autopas::Newton3Option::disabled,
                                                      autopas::InteractionTypeOption::pairwise,
                                                      autopas::VectorizationPatternOption::NA};

  // Triwise Verlet Lists Config
  const autopas::Configuration _confVl_list_iteration_3b{autopas::ContainerOption::verletLists,
                                                         _cellSizeFactor,
                                                         autopas::TraversalOption::vl_list_iteration,
                                                         autopas::LoadEstimatorOption::none,
                                                         autopas::DataLayoutOption::aos,
                                                         autopas::Newton3Option::disabled,
                                                         autopas::InteractionTypeOption::triwise,
                                                         autopas::VectorizationPatternOption::NA};

  // Triwise Verlet Lists Config with Pair List Traversal
  const autopas::Configuration _confVl_pair_list_iteration_3b{autopas::ContainerOption::verletLists,
                                                              _cellSizeFactor,
                                                              autopas::TraversalOption::vl_pair_list_iteration,
                                                              autopas::LoadEstimatorOption::none,
                                                              autopas::DataLayoutOption::aos,
                                                              autopas::Newton3Option::disabled,
                                                              autopas::InteractionTypeOption::triwise,
                                                              autopas::VectorizationPatternOption::NA};
};
