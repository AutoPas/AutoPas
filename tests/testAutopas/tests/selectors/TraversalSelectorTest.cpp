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
TEST_F(TraversalSelectorTest, testSelectAndGetCurrentTraversal) {
  MFunctor functor;

  // this should be high enough so that sliced is still valid for the current processors thread count.
  constexpr size_t domainSize = 1000;
  autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize});

  // expect an exception if nothing is selected yet
  EXPECT_THROW((traversalSelector.generateTraversal<MFunctor, false, false>(autopas::TraversalOption(-1), functor)),
               autopas::utils::ExceptionHandler::AutoPasException);

  for (auto &traversalOption : autopas::allTraversalOptions) {
    auto traversal = traversalSelector.generateTraversal<MFunctor, false, false>(traversalOption, functor);

    // check that traversals are of the expected type
    EXPECT_EQ(traversalOption, traversal->getTraversalType())
        << "Is the domain size large enough for the processors' thread count?";
  }
}

// TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestAbs) {
//  auto strategy = autopas::SelectorStrategy::fastestAbs;
//
//  mapOptionsTime measurements;
//
//  measurements[autopas::TraversalOption::c08] = {22, 14};
//  measurements[autopas::TraversalOption::sliced] = {30, 10};
//
//  mapOptionsTime ignoredMeasurements;
//  ignoredMeasurements[autopas::TraversalOption::c08] = {1};
//
//  testFastest(strategy, measurements, autopas::TraversalOption::sliced, ignoredMeasurements);
//}
//
// TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestMean) {
//  auto strategy = autopas::SelectorStrategy::fastestMean;
//
//  mapOptionsTime measurements;
//
//  measurements[autopas::TraversalOption::c08] = {2, 20};
//  measurements[autopas::TraversalOption::sliced] = {5, 7};
//
//  testFastest(strategy, measurements, autopas::TraversalOption::sliced);
//}
//
// TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestMedian) {
//  auto strategy = autopas::SelectorStrategy::fastestMedian;
//
//  mapOptionsTime measurements;
//
//  measurements[autopas::TraversalOption::c08] = {4, 1, 5};
//  measurements[autopas::TraversalOption::sliced] = {2, 3, 3, 100};
//
//  testFastest(strategy, measurements, autopas::TraversalOption::sliced);
//}
//
// void TraversalSelectorTest::testFastest(autopas::SelectorStrategy strategy, mapOptionsTime measurements,
//                                        autopas::TraversalOption expectedBest, mapOptionsTime ignoredMeasurements) {
//  MFunctor functor;
//
//  std::vector<autopas::TraversalOption> optionVector;
//  optionVector.reserve(measurements.size());
//
//  for (auto &&m : measurements) {
//    optionVector.push_back(m.first);
//  }
//
//  constexpr size_t domainSize = 1000;
//  autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize}, optionVector);
//
//  EXPECT_THROW((traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor)), std::exception);
//
//  for (auto &&m : measurements) {
//    for (auto &&t : m.second) {
//      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(true));
//      traversalSelector.addTimeMeasurement(functor, m.first, t);
//    }
//  }
//
//  for (auto &&m : ignoredMeasurements) {
//    for (auto &&t : m.second) {
//      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(false));
//      traversalSelector.addTimeMeasurement(functor, m.first, t);
//    }
//  }
//
//  auto traversal = traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor);
//  EXPECT_EQ(expectedBest, traversal->getTraversalType());
//
//  // select optimal traversal should delete all measurements
//  EXPECT_THROW((traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor)), std::exception);
//}
