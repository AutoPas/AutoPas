/**
 * @file FeatureVectorTest.h
 * @author Jan Nguyen
 * @date 11.08.19
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/tuning/utils/FeatureVectorEncoder.h"

class FeatureVectorTest : public AutoPasTestBase {
 public:
  FeatureVectorTest();

  std::vector<autopas::FeatureVector::ContainerTraversalEstimatorOption> allCompatibleContainerTraversalEstimators;
  std::vector<autopas::DataLayoutOption> allDataLayouts;
  std::vector<autopas::Newton3Option> allNewton3;
};
