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
  const autopas::Configuration configurationC01 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  const autopas::Configuration configurationC08 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  const autopas::Configuration configurationSliced = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  static constexpr double relativeOptimumRange{1.2};
  static constexpr unsigned int maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int evidenceFirstPrediction{2};
  static constexpr autopas::ExtrapolationMethodOption linePrediction{
      autopas::ExtrapolationMethodOption::linePrediction};
};
