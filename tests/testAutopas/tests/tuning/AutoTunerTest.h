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

  const autopas::Configuration _confLc_c01{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                           autopas::TraversalOption::lc_c01,      autopas::LoadEstimatorOption::none,
                                           autopas::DataLayoutOption::aos,        autopas::Newton3Option::disabled};
  const autopas::Configuration _confLc_c18{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                           autopas::TraversalOption::lc_c18,      autopas::LoadEstimatorOption::none,
                                           autopas::DataLayoutOption::aos,        autopas::Newton3Option::disabled};
  const autopas::Configuration _confLc_c08{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                           autopas::TraversalOption::lc_c08,      autopas::LoadEstimatorOption::none,
                                           autopas::DataLayoutOption::aos,        autopas::Newton3Option::disabled};
};
