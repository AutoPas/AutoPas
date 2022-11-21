/**
 * @file AutoTunerTest.h
 * @author F. Gratl
 * @date 8/10/18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/selectors/Configuration.h"

class AutoTunerTest : public AutoPasTestBase {
 public:
  AutoTunerTest() = default;
  ~AutoTunerTest() override = default;

  const double _cellSizeFactor{1.};

  const autopas::Configuration _confLc_c01{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                           autopas::TraversalOption::lc_c01,      autopas::LoadEstimatorOption::none,
                                           autopas::DataLayoutOption::aos,        autopas::Newton3Option::disabled, 5};
  const autopas::Configuration _confLc_c04{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                           autopas::TraversalOption::lc_c04,      autopas::LoadEstimatorOption::none,
                                           autopas::DataLayoutOption::aos,        autopas::Newton3Option::disabled, 5};
  const autopas::Configuration _confLc_c08{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                           autopas::TraversalOption::lc_c08,      autopas::LoadEstimatorOption::none,
                                           autopas::DataLayoutOption::aos,        autopas::Newton3Option::disabled, 5};
};
