/**
 * @file TuningTest.cpp
 * @author S. Newcome
 * @date 29.10.2026
 */

#include "TuningTest.h"

#include <optional>
#include <set>

#include "autopas/AutoPasDecl.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(MockPairwiseFunctor<Molecule> *);

/**
 * Tests that a tuning phase which ends on a rejected configuration behaves correctly.
 *
 * The search space consists of lc_c08 and lc_c18 with and without Newton3, with a mock functor that does not support
 * Newton3, so both Newton3 configurations are rejected. The last configuration tested in the first tuning phase is one
 * of them.
 */
TEST_F(TuningTest, testEndingTuningPhaseWithRejectedConfig) {
  constexpr size_t numSamples = 2;
  constexpr size_t tuningInterval = 20;
  // Only the two configurations without Newton3 are applicable.
  constexpr size_t numApplicableConfigs = 2;
  constexpr size_t numTuningIterations = numApplicableConfigs * numSamples;

  autopas::AutoPas<Molecule> autoPas;
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
  autoPas.setAllowedTraversals({autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c18});
  autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled, autopas::Newton3Option::enabled});
  autoPas.setNumSamples(numSamples);
  autoPas.setTuningInterval(tuningInterval);
  autoPas.init();
  ASSERT_FALSE(autoPas.searchSpaceIsTrivial());

  size_t iteration = 0;
  // Iterations in which at least one configuration was rejected.
  std::set<size_t> iterationsWithRejection{};

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  // AutoPas asks the functor before using a Newton3 configuration if allowed, so every call is a rejection.
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly([&]() {
    iterationsWithRejection.insert(iteration);
    return false;
  });
  MockPairwiseFunctor<Molecule> *functorPtr = &functor;

  // Runs one tuning interval and checks that only the first numTuningIterations iterations are tuning iterations, that
  // no rejected configuration is ever used, and that the selected configuration does not change after tuning.
  const auto runTuningInterval = [&]() {
    const size_t firstIteration = iteration;
    std::optional<autopas::Configuration> selectedConfig{};
    for (; iteration < firstIteration + tuningInterval; ++iteration) {
      [[maybe_unused]] const auto leavingParticles = autoPas.updateContainer();
      const bool stillTuning = autoPas.computeInteractions(functorPtr);
      const autopas::Configuration config = autoPas.getCurrentConfigs().at(autopas::InteractionTypeOption::pairwise);

      const bool expectTuning = iteration < firstIteration + numTuningIterations;
      EXPECT_EQ(stillTuning, expectTuning) << "Iteration " << iteration;
      EXPECT_EQ(config.newton3, autopas::Newton3Option::disabled)
          << "Iteration " << iteration << " used the rejected configuration " << config;

      // If not tuning, the config should not change
      if (not expectTuning) {
        if (not selectedConfig.has_value()) {
          selectedConfig = config;
        }
        EXPECT_EQ(config, selectedConfig.value())
            << "The selected configuration changed outside of the tuning phase in iteration " << iteration;
      }
    }
  };

  // First tuning phase, which ends on a rejected configuration.
  runTuningInterval();
  // Make sure that the first tuning phase really ended on a rejected configuration, i.e. that a configuration was
  // rejected in the first iteration after the tuning phase.
  ASSERT_TRUE(iterationsWithRejection.contains(numTuningIterations))
      << "The last configuration of the first tuning phase was not rejected, so this test does not test what it "
         "should.";

  // Second tuning phase. Tuning has to continue as usual after a tuning phase that ended on a rejected configuration.
  runTuningInterval();
}