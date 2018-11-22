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

TEST_F(TraversalSelectorTest, testSelectOptimalTraversal) {
  MFunctor functor;

  std::vector<autopas::TraversalOptions> optionVector = {autopas::TraversalOptions::sliced,
                                                         autopas::TraversalOptions::c08};

  constexpr size_t domainSize = 1000;
  autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize}, optionVector);

  EXPECT_THROW(
      (traversalSelector.selectOptimalTraversal<MFunctor, true, true>(autopas::SelectorStrategy::fastestAbs, functor)),
      std::exception);

  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(Return(true));
  traversalSelector.addTimeMeasurement(functor, optionVector[0], 20);
  traversalSelector.addTimeMeasurement(functor, optionVector[0], 22);
  traversalSelector.addTimeMeasurement(functor, optionVector[1], 30);
  traversalSelector.addTimeMeasurement(functor, optionVector[1], 10);

  // add one really fast time which is ignored
  EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(false));
  traversalSelector.addTimeMeasurement(functor, optionVector[0], 1);

  auto traversal =
      traversalSelector.selectOptimalTraversal<MFunctor, true, true>(autopas::SelectorStrategy::fastestAbs, functor);
  EXPECT_EQ(optionVector[1], traversal->getTraversalType());

  // select optimal traversal should delete all measurements
  EXPECT_THROW(
      (traversalSelector.selectOptimalTraversal<MFunctor, true, true>(autopas::SelectorStrategy::fastestAbs, functor)),
      std::exception);
}