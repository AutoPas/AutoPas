/**
 * @file AutoPasConfigurationCommunicatorTest.h
 * @author muehlhaeusser
 * @date 16.05.24
 */

#pragma once

#include <gtest/gtest.h>

#include <set>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/utils/NumberSetFinite.h"

class AutoPasConfigurationCommunicatorTest : public AutoPasTestBase {
 public:
  void testConfigsCommunication(std::set<autopas::Configuration> &configs);

  std::set<autopas::ContainerOption> containerOptions = autopas::ContainerOption::getMostOptions();
  std::set<autopas::LoadEstimatorOption> loadEstimatorOptions = autopas::LoadEstimatorOption::getMostOptions();
  std::set<autopas::DataLayoutOption> dataLayoutOptions = autopas::DataLayoutOption::getMostOptions();
  std::set<autopas::Newton3Option> newton3Options = autopas::Newton3Option::getMostOptions();
  autopas::NumberSetFinite<double> cellSizeFactors = {0.5, 1.0, 1.5};
  std::set<autopas::TraversalOption> pairwiseTraversalOptions = autopas::TraversalOption::getMostPairwiseOptions();
  std::set<autopas::TraversalOption> triwiseTraversalOptions = autopas::TraversalOption::getMostTriwiseOptions();
};