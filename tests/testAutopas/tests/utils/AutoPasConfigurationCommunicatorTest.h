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
  void testConfigsCommunication(const std::set<autopas::Configuration> &configs);

  const std::set<autopas::ContainerOption> containerOptions = autopas::ContainerOption::getMostOptions();
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions = autopas::LoadEstimatorOption::getMostOptions();
  const std::set<autopas::DataLayoutOption> dataLayoutOptions = autopas::DataLayoutOption::getMostOptions();
  const std::set<autopas::Newton3Option> newton3Options = autopas::Newton3Option::getMostOptions();
  const autopas::NumberSetFinite<double> cellSizeFactors = {0.5, 1.0, 1.5};
  const std::set<autopas::TraversalOption> pairwiseTraversalOptions =
      autopas::TraversalOption::getMostPairwiseOptions();
  const std::set<autopas::TraversalOption> triwiseTraversalOptions = autopas::TraversalOption::getMostTriwiseOptions();
  const std::set<autopas::VectorizationPatternOption> vecPatternOptions =
      autopas::VectorizationPatternOption::getMostOptions();
};