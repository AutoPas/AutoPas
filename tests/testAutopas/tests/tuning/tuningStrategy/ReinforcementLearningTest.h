/**
 * @file ReinforcementLearningTest.h
 * @author p. Metscher
 * @date 24.06.25
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/tuning/tuningStrategy/ReinforcementLearning.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"

class ReinforcementLearningTest : public AutoPasTestBase {
 private:
  /**
   * The container options for the cartesian product search space.
   */
  const std::set<autopas::ContainerOption> containerOptions = autopas::ContainerOption::getMostOptions();

  /**
   * The traversal options for the cartesian product search space.
   */
  const std::set<autopas::TraversalOption> traversalOptions = autopas::TraversalOption::getMostOptions();

  /**
   * The load estimator options for the cartesian product search space.
   */
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions = autopas::LoadEstimatorOption::getMostOptions();

  /**
   * The data layout options for the cartesian product search space.
   */
  const std::set<autopas::DataLayoutOption> dataLayoutOptions = autopas::DataLayoutOption::getMostOptions();

  /**
   * The newton3 options for the cartesian product search space.
   */
  const std::set<autopas::Newton3Option> newton3Options = autopas::Newton3Option::getMostOptions();

  /**
   * The cell size factors for the cartesian product search space.
   */
  const autopas::NumberSetFinite<double> cellSizeFactors{0.5, 1.0};

 public:
  /**
   * The full cartesian product serving as the base configuration queue.
   */
  const std::set<autopas::Configuration> configQueueBase = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      autopas::InteractionTypeOption::pairwise);
};
