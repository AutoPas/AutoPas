/**
 * @file PredictiveTuningTest.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/tuningStrategy/PredictiveTuning.h"

class PredictiveTuningTest : public AutoPasTestBase,
                             public ::testing::WithParamInterface<autopas::ExtrapolationMethodOption> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto extrapolationOption = static_cast<ParamType>(info.param);
      auto str = extrapolationOption.to_string();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };

 protected:
  /**
   * Tunes for a few iterations (length of evidence) and checks whether all possible configurations were tuned.
   * @param predictiveTuning The PredictiveTuning strategy.
   * @param evidence Map of configurations to evidence. These simulate the performance of every configuration.
   * @param expectedPredictions Map of configurations to expected predictions. Pass an empty map if this should be
   * skipped.
   * @param iterationCounter The current iteration number; will be increased accordingly.
   * @param tuningPhase Number of the current tuningPhase.
   * @return Vector of the tested configurations.
   */
  static auto simulateTuningPhase(autopas::PredictiveTuning &predictiveTuning,
                                  const std::map<autopas::Configuration, long> &newEvidence,
                                  autopas::EvidenceCollection &evidenceCollection, size_t &iterationCounter,
                                  size_t tuningPhase) {
    // input sanity checks
    ASSERT_FALSE(newEvidence.empty()) << "There should be predetermined evidence.";

    // fill the queue with everthing that is in newEvidence
    std::vector<autopas::Configuration> confQueue;
    for (const auto &[conf, _] : newEvidence) {
      confQueue.push_back(conf);
    }
    const auto configsExpectedToBeTested = confQueue;
    // collects all configurations that are tested in this phase
    std::vector<autopas::Configuration> testedConfigs;
    testedConfigs.reserve(confQueue.size());

    // strategies are reset at the start of every phase.
    predictiveTuning.reset(iterationCounter, tuningPhase, confQueue, evidenceCollection);

    // identify optimum of the previous tuning phase and check that this is now at the back of the queue
    if (tuningPhase > _evidenceFirstPrediction) {
      const auto [bestConfigLastPhase, _] = evidenceCollection.getOptimalConfiguration(tuningPhase - 1);
      EXPECT_EQ(confQueue.back(), bestConfigLastPhase)
          << "We expect the best predicted configuration to be at the back of the "
             "FiLo queue at the start of the tuning phase. Tuning Phase: "
          << tuningPhase;
    }

    // simulate the iterations of the tuning phase
    while (not confQueue.empty()) {
      const auto &conf = confQueue.back();
      testedConfigs.push_back(conf);

      // find this configuration's evidence in the planned evidence
      EXPECT_NO_THROW(newEvidence.at(conf)) << "The newEvidence has no data for " << conf << ".\nCheck the test setup!";

      const autopas::Evidence evidence{iterationCounter, tuningPhase, newEvidence.at(conf)};
      predictiveTuning.addEvidence(conf, evidence);
      evidenceCollection.addEvidence(conf, evidence);
      confQueue.pop_back();  // this line invalidates conf
      predictiveTuning.optimizeSuggestions(confQueue, evidenceCollection);
      ++iterationCounter;
    }

    EXPECT_THAT(testedConfigs, testing::UnorderedElementsAreArray(configsExpectedToBeTested))
        << "Test did not test all expected configurations in tuning phase.";

    // at the end of a tuning phase the strategy should point to the optimum
    const auto [bestConfigThisPhase, _] =
        *(std::min_element(newEvidence.begin(), newEvidence.end(),
                           [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; }));
    EXPECT_EQ(std::get<0>(evidenceCollection.getLatestOptimalConfiguration()), bestConfigThisPhase)
        << "Did not select the fastest configuration. At the end of the tuning phase.";
  }

  /**
   * TODO: Update doc
   * Simulates multiple tuning phases and checks if the predictions, inferred from the given evidence matches the
   * expectations.
   * The number of phases depends on the length of the provided vector of evidence.
   *
   * A predictive tuning object is generated and a configQueue, which is filled with all configurations from the
   * expectations map.
   * The provided evidence is fed to the tuning strategy as if it was tuning for several phases.
   *
   * @param extrapolationMethodOption
   * @param evidencePerPhase
   * @param tuningInterval Number of iterations between the last tuning iteration and evaluation of expectations.
   * @param expectedPredictions
   */
  static void checkPredictions(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                               const std::vector<std::map<autopas::Configuration, long>> &evidencePerPhase,
                               unsigned int tuningInterval,
                               const std::map<autopas::Configuration, long> &expectedPredictions);

  static constexpr autopas::Configuration _configurationLC_C01 =
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
                             autopas::Configuration::ThreadCountNoTuning);
  static constexpr autopas::Configuration _configurationLC_C08 =
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
                             autopas::Configuration::ThreadCountNoTuning);

  static constexpr autopas::Configuration _configurationLC_Sliced =
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
                             autopas::Configuration::ThreadCountNoTuning);

  static constexpr double _relativeOptimumRange{1.2};
  static constexpr unsigned int _maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int _evidenceFirstPrediction{2};
};
