/**
 * @file PredictiveTuningTest.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/selectors/tuningStrategy/PredictiveTuning.h"

class PredictiveTuningTest : public AutoPasTestBase {
 protected:
  const autopas::Configuration configurationLC_C01 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  const autopas::Configuration configurationLC_C08 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  const autopas::Configuration configurationLC_Sliced = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
      autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  std::vector<autopas::Configuration> allConfigs{configurationLC_C01, configurationLC_C08, configurationLC_Sliced};

  static constexpr double relativeOptimumRange{1.2};
  static constexpr unsigned int maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int evidenceFirstPrediction{2};
  static constexpr autopas::ExtrapolationMethodOption linePrediction{
      autopas::ExtrapolationMethodOption::linePrediction};
};
