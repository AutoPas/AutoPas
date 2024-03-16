/**
 * @file PredictiveTuningTest.cpp
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#include "PredictiveTuningTest.h"

#include <gmock/gmock-more-matchers.h>

#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/PredictiveTuning.h"

void PredictiveTuningTest::checkPredictions(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                                            const std::vector<std::map<autopas::Configuration, long>> &evidencePerPhase,
                                            unsigned int tuningInterval,
                                            const std::map<autopas::Configuration, long> &expectedPredictions) {
  // setup sanity check
  for (const auto &phaseEvidence : evidencePerPhase) {
    ASSERT_EQ(phaseEvidence.size(), expectedPredictions.size())
        << "Provided number of evidence and expectation maps do not match! "
           "The test needs one prediction for each configuration that gets evidences until the last phase.";
  }

  // infer the configurations to be tested from the expected evidence
  std::vector<autopas::Configuration> confQueue{};
  confQueue.reserve(expectedPredictions.size());
  for (const auto &[conf, _] : expectedPredictions) {
    confQueue.push_back(conf);
  }
  const auto configsExpectedToBeTested = confQueue;

  autopas::PredictiveTuning predictiveTuning{_relativeOptimumRange, _maxTuningIterationsWithoutTest,
                                             static_cast<unsigned int>(evidencePerPhase.size()),
                                             extrapolationMethodOption};

  // insert evidence into collection and count in which iteration and tuning phase we are
  size_t iteration = 0;
  size_t tuningPhase = 0;
  autopas::EvidenceCollection evidenceCollection{};
  for (const auto &evidenceMap : evidencePerPhase) {
    for (const auto &[conf, value] : evidenceMap) {
      const autopas::Evidence evidence{iteration, tuningPhase, value};
      evidenceCollection.addEvidence(conf, evidence);
      predictiveTuning.addEvidence(conf, evidence);
      ++iteration;
    }
    iteration += tuningInterval;
    ++tuningPhase;
  }

  // calculate predictions for the subsequent tuning phase
  const auto predictions =
      predictiveTuning.calculatePredictions(iteration, tuningPhase, configsExpectedToBeTested, evidenceCollection);

  // compare predictions to expectations.
  for (const auto &[conf, pred] : predictions) {
    EXPECT_EQ(pred, expectedPredictions.at(conf))
        << "Expected and calculated predictions differ for " << conf.toString();
  }
}

/**
 * Tests the prediction and the selection of the right optimum.
 * Three different configurations:
 *      - first that starts optimal but becomes slowly more expensive.
 *      - second starts expensive but becomes quickly better. (expected optimum)
 *      - third is constant.
 * @note All expected values in this test were calculated by hand and are precise if not stated otherwise.
 */
TEST_P(PredictiveTuningTest, testPredictions) {
  auto extrapolationOption = GetParam();
  // we do the generation + switch case here to ensure a test is generated for every method.
  switch (extrapolationOption) {
    case autopas::ExtrapolationMethodOption::linePrediction: {
      checkPredictions(extrapolationOption,
                       {
                           {{_configurationLC_C01, 94}, {_configurationLC_C08, 109}, {_configurationLC_Sliced, 101}},
                           {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}},
                       },
                       // all predictions are evaluated for the seventh iteration (iteration==6)
                       0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 99}, {_configurationLC_Sliced, 101}});
      break;
    }
    case autopas::ExtrapolationMethodOption::newton: {
      checkPredictions(extrapolationOption,
                       {
                           {{_configurationLC_C01, 79}, {_configurationLC_C08, 115}, {_configurationLC_Sliced, 101}},
                           {{_configurationLC_C01, 90}, {_configurationLC_C08, 109}, {_configurationLC_Sliced, 101}},
                           {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}},
                       },
                       // all predictions are evaluated for the tenth iteration (iteration==9)
                       // for _configurationLC_C08 we actually expect 99.33
                       0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 99}, {_configurationLC_Sliced, 101}});
      break;
    }
    case autopas::ExtrapolationMethodOption::linearRegression: {
      checkPredictions(
          extrapolationOption,
          {
              // values along the desired line offset so that the regression should still return the same line
              {{_configurationLC_C01, 91 + 1}, {_configurationLC_C08, 115 - 1}, {_configurationLC_Sliced, 101 - 2}},
              {{_configurationLC_C01, 94 - 2}, {_configurationLC_C08, 109 + 2}, {_configurationLC_Sliced, 101 + 4}},
              {{_configurationLC_C01, 97 + 1}, {_configurationLC_C08, 103 - 1}, {_configurationLC_Sliced, 101 - 2}},
          },
          // all predictions are evaluated for the tenth iteration (iteration==9)
          0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 99}, {_configurationLC_Sliced, 101}});
      break;
    }
    case autopas::ExtrapolationMethodOption::lastResult: {
      checkPredictions(extrapolationOption,
                       {
                           {{_configurationLC_C01, 94}, {_configurationLC_C08, 109}, {_configurationLC_Sliced, 101}},
                           {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}},
                       },
                       // all predictions are evaluated for the seventh iteration (iteration==6)
                       0, {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}});
      break;
    }
    default:
      GTEST_FAIL() << "No test implemented for extrapolation method " << extrapolationOption.to_string();
  }
}

/**
 * Tests the prediction and the selection of the right optimum. Values are chosen so that the predictions result in
 * under or overflow.
 *
 * @note We expect any underflow to be set to 1 and overflows to max size_t - 1.
 *
 * (max-value itself is used as a indicator that no prediction was computed).
 */
TEST_P(PredictiveTuningTest, testUnderAndOverflow) {
  auto extrapolationOption = GetParam();
  constexpr auto maxLong = std::numeric_limits<long>::max();
  std::map<autopas::Configuration, long> expected = {
      {_configurationLC_C01, 1}, {_configurationLC_C08, maxLong - 1}, {_configurationLC_Sliced, 101}};
  if (extrapolationOption == autopas::ExtrapolationMethodOption::lastResult) {
    expected = {{_configurationLC_C01, 10}, {_configurationLC_C08, maxLong - 100}, {_configurationLC_Sliced, 101}};
  }
  checkPredictions(
      extrapolationOption,
      {
          {{_configurationLC_C01, 100}, {_configurationLC_C08, 1}, {_configurationLC_Sliced, 101}},
          {{_configurationLC_C01, 10}, {_configurationLC_C08, maxLong - 100}, {_configurationLC_Sliced, 101}},
      },
      // all predictions are evaluated for the seventh iteration (iteration==6)
      100, expected);
}

INSTANTIATE_TEST_SUITE_P(Generated, PredictiveTuningTest,
                         ::testing::ValuesIn(autopas::ExtrapolationMethodOption::getAllOptions()),
                         PredictiveTuningTest::PrintToStringParamName());

/**
 * Tests the selection of the right optimumSearchSpace in the regard of the number of iterations without a test and
 * tuning. Two different configurations:
 *      - first is constantly outside of the optimum range (20).
 *      - second is constantly the optimum (10).
 * Tuning phases 1 and 2 should cover the whole search space bc we need at least two samples to make predictions
 * In tuning phase 3-8 only the second configuration should be in _optimalSearchSpace.
 * In tuning phase 9  both configuration should be tested because the slow should be reinserted via the
 * tooLongNotTested mechanic.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTooLongNotTested) {
  size_t iteration = 0;
  size_t tuningPhase = 0;
  const std::set<autopas::Configuration> searchSpace{_configurationLC_C08, _configurationLC_Sliced};
  autopas::EvidenceCollection evidenceCollection{};

  auto predictiveTuning =
      autopas::PredictiveTuning(_relativeOptimumRange, _maxTuningIterationsWithoutTest, _evidenceFirstPrediction,
                                autopas::ExtrapolationMethodOption::linePrediction);

  std::map<autopas::Configuration, long> evidenceBest{{_configurationLC_C08, 10}};
  // Not a hard requirement but if this is violated the merge into evidenceAll needs to be done differently.
  ASSERT_EQ(evidenceBest.size(), 1);
  std::map<autopas::Configuration, long> evidenceAll{{_configurationLC_C01, 20}, *evidenceBest.begin()};

  std::vector<std::pair<autopas::Configuration, long>> evidenceSorted(evidenceAll.begin(), evidenceAll.end());
  std::sort(evidenceSorted.begin(), evidenceSorted.end(),
            [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; });

  // PHASE A: Initial tuning phases where every configuration is tested because we can not make predictions yet
  for (; tuningPhase < _evidenceFirstPrediction; ++tuningPhase) {
    simulateTuningPhase(predictiveTuning, evidenceAll, evidenceCollection, iteration, tuningPhase);
  }
  // Check for expected number of iterations
  const auto expectedNumIterationsInitial = _evidenceFirstPrediction * evidenceAll.size();
  EXPECT_EQ(iteration, expectedNumIterationsInitial) << "First tuning phases did not test all configurations!";

  // PHASE B: Predictions are made and non-near-optimal configurations are not tested
  for (; tuningPhase < _evidenceFirstPrediction + _maxTuningIterationsWithoutTest; ++tuningPhase) {
    // We pass only evidence information for configurations that should be tested. If the tuning strategy wants to
    // test anything else we get an exception that the requested conf is not in the evidence map.
    EXPECT_NO_THROW(simulateTuningPhase(predictiveTuning, evidenceBest, evidenceCollection, iteration, tuningPhase))
        << "Iteration: " << iteration << "\nTesting a configuration other than " << evidenceBest.begin()->first;
  }
  // Check for expected number of iterations
  const auto expectedNumIterationsPredicting = _maxTuningIterationsWithoutTest * evidenceBest.size();
  EXPECT_EQ(iteration, expectedNumIterationsInitial + expectedNumIterationsPredicting)
      << "Tuning phases with predictions tested unexpected number of configurations!";

  // PHASE C: Now the non-near-optimal configurations should be tested again
  simulateTuningPhase(predictiveTuning, evidenceAll, evidenceCollection, iteration, tuningPhase);
  // Check for expected number of iterations
  EXPECT_EQ(iteration, expectedNumIterationsInitial + expectedNumIterationsPredicting + evidenceAll.size())
      << "Non-Optimal Conf is not retested after cool down time.";
}
