/**
 * @file C01TraversalAdaptiveTest.cpp
 * @author C. Menges
 * @date 11.07.2019
 */

#include "C01TraversalAdaptiveTest.h"

TEST_F(C01TraversalAdaptiveTest, testIsApplicable) {
  MFunctor mFunctor;
  {
    std::unique_ptr<autopas::TraversalSelectorInfoAdaptive<FPCell>> info;
    autopas::C01TraversalAdaptive<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> C01TraversalAdaptiveN3On(
        &mFunctor, std::move(info));
    EXPECT_EQ(C01TraversalAdaptiveN3On.isApplicable(), false);
  }
  {
    std::unique_ptr<autopas::TraversalSelectorInfoAdaptive<FPCell>> info;
    autopas::C01TraversalAdaptive<FPCell, MFunctor, autopas::DataLayoutOption::aos, false> C01TraversalAdaptiveN3On(
        &mFunctor, std::move(info));
    EXPECT_EQ(C01TraversalAdaptiveN3On.isApplicable(), true);
  }
}
