/**
 * @file TraversalSelectorTest.cpp
 * @author F. Gratl
 * @date 21.06.18
 */

#include "TraversalSelectorTest.h"

using ::testing::Return;

/**
 * Check if the only allowed traversal is returned
 */
TEST_F(TraversalSelectorTest, testGetOptimalTraversalOneOption) {
  MFunctor functor;

  for (auto &traversalOption : autopas::allTraversalOptions) {
    // this should be high enough so that sliced is still valid for the current processors thread count.
    constexpr size_t domainSize = 1000;

    autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize}, {traversalOption});

    EXPECT_THROW((traversalSelector.getOptimalTraversal<MFunctor, false, true>(functor)), std::exception);

    traversalSelector.selectNextTraversal<MFunctor, false, true>(functor);

    auto traversal = traversalSelector.getOptimalTraversal<MFunctor, false, true>(functor);

    // check that traversals are of the expected type
    EXPECT_EQ(traversalOption, traversal->getTraversalType())
        << "Is the domain size large enough for the processors' thread count?";

    // now that the functor is known check if still the same is returned
    traversal = traversalSelector.getOptimalTraversal<MFunctor, false, true>(functor);
    // check that traversals are of the expected type
    EXPECT_EQ(traversalOption, traversal->getTraversalType())
        << "Repeated call for traversal " << traversalOption << " failed";
  }
}

TEST_F(TraversalSelectorTest, testGetOptimalTraversalBadFirstOption) {
  MFunctor functor;

  std::vector<autopas::TraversalOptions> optionVector = {autopas::TraversalOptions::sliced,
                                                         autopas::TraversalOptions::c08};

  autopas::TraversalSelector<FPCell> traversalSelectorC08({1, 1, 1}, optionVector);
  EXPECT_THROW((traversalSelectorC08.getOptimalTraversal<MFunctor, false, true>(functor)), std::exception);
  traversalSelectorC08.selectNextTraversal<MFunctor, false, true>(functor);
  auto traversal = traversalSelectorC08.getOptimalTraversal<MFunctor, false, true>(functor);

  // check that traversals are of the expected type
  EXPECT_EQ(autopas::TraversalOptions::c08, traversal->getTraversalType());

  // also after the functor is known
  traversal = traversalSelectorC08.getOptimalTraversal<MFunctor, false, true>(functor);
  EXPECT_EQ(autopas::TraversalOptions::c08, traversal->getTraversalType());
}

TEST_F(TraversalSelectorTest, testNextTraversal) {
  MFunctor functor;

  std::vector<autopas::TraversalOptions> optionVector = {autopas::TraversalOptions::sliced,
                                                         autopas::TraversalOptions::c08};

  constexpr size_t domainSize = 1000;
  autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize}, optionVector);

  // assure that this works repeatedly
  for (int i = 0; i < 2; ++i) {
    auto errorText = "Failed in repetition " + std::to_string(i);
    auto traversal = traversalSelector.selectNextTraversal<MFunctor, true, true>(functor);
    EXPECT_EQ(autopas::TraversalOptions::sliced, traversal->getTraversalType()) << errorText;

    traversal = traversalSelector.selectNextTraversal<MFunctor, true, true>(functor);
    EXPECT_EQ(autopas::TraversalOptions::c08, traversal->getTraversalType()) << errorText;

    traversal = traversalSelector.selectNextTraversal<MFunctor, true, true>(functor);
    EXPECT_EQ(nullptr, traversal) << errorText;
  }
}

TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestAbs) {
  auto strategy = autopas::SelectorStrategy::fastestAbs;

  std::unordered_map<autopas::TraversalOptions, std::vector<long>> measurements;

  measurements[autopas::TraversalOptions::c08] = {22, 14};
  measurements[autopas::TraversalOptions::sliced] = {30, 10};

  std::unordered_map<autopas::TraversalOptions, std::vector<long>> ignoredMeasurements;
  ignoredMeasurements[autopas::TraversalOptions::c08] = {1};

  testFastest(strategy, measurements, autopas::TraversalOptions::sliced, ignoredMeasurements);
}

TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestMean) {
  auto strategy = autopas::SelectorStrategy::fastestMean;

  std::unordered_map<autopas::TraversalOptions, std::vector<long>> measurements;

  measurements[autopas::TraversalOptions::c08] = {2, 20};
  measurements[autopas::TraversalOptions::sliced] = {5, 7};

  testFastest(strategy, measurements, autopas::TraversalOptions::sliced);
}

TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestMedian) {
  auto strategy = autopas::SelectorStrategy::fastestMedian;

  std::unordered_map<autopas::TraversalOptions, std::vector<long>> measurements;

  measurements[autopas::TraversalOptions::c08] = {4, 1, 5};
  measurements[autopas::TraversalOptions::sliced] = {2, 3, 3, 100};

  testFastest(strategy, measurements, autopas::TraversalOptions::sliced);
}

void TraversalSelectorTest::testFastest(
    autopas::SelectorStrategy strategy, std::unordered_map<autopas::TraversalOptions, std::vector<long>> measurements,
    autopas::TraversalOptions expectedBest,
    std::unordered_map<autopas::TraversalOptions, std::vector<long>> ignoredMeasurements) {
  MFunctor functor;

  std::vector<autopas::TraversalOptions> optionVector;
  optionVector.reserve(measurements.size());

  for (auto &&m : measurements) {
    optionVector.push_back(m.first);
  }

  constexpr size_t domainSize = 1000;
  autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize}, optionVector);

  EXPECT_THROW((traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor)), std::exception);

  for (auto &&m : measurements) {
    for (auto &&t : m.second) {
      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(true));
      traversalSelector.addTimeMeasurement(functor, m.first, t);
    }
  }

  for (auto &&m : ignoredMeasurements) {
    for (auto &&t : m.second) {
      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(false));
      traversalSelector.addTimeMeasurement(functor, m.first, t);
    }
  }

  auto traversal = traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor);
  EXPECT_EQ(expectedBest, traversal->getTraversalType());

  // select optimal traversal should delete all measurements
  EXPECT_THROW((traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor)), std::exception);
}
