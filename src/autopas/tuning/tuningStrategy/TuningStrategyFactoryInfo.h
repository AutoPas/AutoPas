/**
 * @file TuningStrategyFactoryInfo.h
 * @author F. Gratl
 * @date 30.06.23
 */

#pragma once

#include <string>

#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ExplorationMethodOption.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/ReinforcementModelOption.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas {

/**
 * Helper struct encapsulating most information needed to build TuningStrategies by the TuningStrategyFactory.
 * This way if only a specific option should be built, not all arguments have to be specified explicitly.
 */
struct TuningStrategyFactoryInfo {
  // Used by multiple Strategies
  /**
   * Strategies that don't converge (or not quickly enough) can be told to limit the number of evidence to this number.
   */
  unsigned int maxEvidence{10};
  // Predictive Tuning Options
  /**
   * Function option used for extrapolating performance from observed evidence.
   */
  ExtrapolationMethodOption extrapolationMethodOption{ExtrapolationMethodOption::linearRegression};
  /**
   * Factor of the range of the optimal configurations for the optimalSearchSpace.
   */
  double relativeOptimumRange{1.2};
  /**
   * If a config is not tested for this number of tuning phases test it again to make predictions more reliable.
   */
  unsigned int maxTuningPhasesWithoutTest{100};
  /**
   * The number of evidence that have to be collected until the first prediction can be made.
   */
  unsigned int minNumberOfEvidence{3};

  // Slow Config Filter Options
  /**
   * Any configuration that is slower than the fastest times this factor will be blacklisted.
   */
  double relativeBlacklistRange{3};

  // Bayesian Strategies Options
  /**
   * Function used to predict informational gain.
   */
  AcquisitionFunctionOption acquisitionFunctionOption{AcquisitionFunctionOption::upperConfidenceBound};

  // Deep Reinforcement Learning Options
  /**
   * The method used for selecting the exploration samples in the deep reinforcement learning tuning strategy.
   */
  ExplorationMethodOption drlExplorationMethod{ExplorationMethodOption::polynomial};
  /**
   * The number of configurations to explore during the exploration phase.
   */
  unsigned int drlRetrainingIterations{34};
  /**
   * The learning rate for the neural network update in the deep reinforcement learning tuning strategy.
   */
  double drlLearningRate{1.771665733218049020e-09};
  /**
   * The number of configurations to explore during the exploration phase.
   */
  unsigned int drlNumExplorationSamples{4};
  /**
   * The number of configurations to explore during the exploitation phase.
   */
  unsigned int drlNumExploitationSamples{1};
  /**
   * The scaling factor used for the priority of the age of a data in polynomial exploration in the deep reinforcement
   * learning tuning strategy.
   */
  double drlPhaseScale{1.162197398643888948e-05};
  /**
   * The weight used for updating the neural network in the reinforcement learning process in the deep reinforcement
   * learning tuning strategy.
   */
  double drlUpdateWeight{0.5};

  // Reinforcement Learning Options
  /**
   * Learning rate for the reinforcement learning algorithm.
   */
  double rlLearningRate{0.8};
  /**
   * Discount factor for the reinforcement learning algorithm.
   */
  double rlDiscountFactor{0.8};
  /**
   * Number of configurations to explore during the exploration phase.
   */
  unsigned int rlExplorationSamples{4};
  /**
   * The model used for the reinforcement learning algorithm.
   */
  ReinforcementModelOption rlModel{ReinforcementModelOption::Value::SARSA};

  // Rule Based Tuning Options
  /**
   * The name and path of the file where the rules are stored.
   */
  std::string ruleFileName{"tuningRules.rule"};

  // Fuzzy Tuning Options
  /**
   * The name and path of the file where the rules are stored.
   */
  std::string fuzzyRuleFileName{"fuzzyRulesSuitability.frule"};

  // MPI Tuning Options
  /**
   * If MPIParallelizedStrategy is in the list of strategies this should be set to true to notify other strategies
   * if necessary.
   */
  bool mpiDivideAndConquer{false};
  /**
   * Maximum absolute difference in similarity metric for two ranks to fall in the same bucket.
   */
  double mpiTuningMaxDifferenceForBucket{0.3};
  /**
   * Weight for maxDensity in the calculation for bucket distribution.
   */
  double mpiTuningWeightForMaxDensity{0.0};
  /**
   * MPI Communicator used within AutoPas.
   */
  AutoPas_MPI_Comm autopasMpiCommunicator{AUTOPAS_MPI_COMM_NULL};
};
}  // namespace autopas
