/**
 * @file TraversalSelectorTest.cpp
 * @author F. Gratl
 * @date 21.06.18
 */

#include "TraversalSelectorTest.h"

#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::Return;

/**
 * Check if the only allowed traversal is returned
 */
TEST_F(TraversalSelectorTest, testSelectAndGetCurrentTraversal) {
  // this should be high enough so that sliced is still valid for the current processors thread count.
  constexpr size_t domainSize = 900;
  autopas::TraversalSelectorInfo traversalSelectorInfo({domainSize, domainSize, domainSize}, 1., {1., 1., 1.}, 32);

  auto checkTraversal = [](const auto &traversalOption, const auto traversal) {
    // check that traversals are of the expected type
    EXPECT_EQ(traversalOption, traversal ? traversal->getTraversalType() : autopas::TraversalOption())
        << "Is the domain size large enough for the processors' thread count?";
  };

  MPairwiseFunctor pairwiseFunctor;
  for (const auto &traversalOption : autopas::TraversalOption::getAllPairwiseOptions()) {
    autopas::DataLayoutOption datalayout = autopas::DataLayoutOption::aos;
    bool newton3 = false;
    if (traversalOption == autopas::TraversalOption::lc_c01_combined_SoA or
        traversalOption == autopas::TraversalOption::lc_c04_combined_SoA) {
      datalayout = autopas::DataLayoutOption::soa;  // these traversals only work with SoA
    }
    if (traversalOption == autopas::TraversalOption::ot_c18) {
      newton3 = true;  // this traversal only supports Newton3 enabled
    }
    auto traversal = autopas::TraversalSelector::generateTraversal<FPCell, MPairwiseFunctor>(
        traversalOption, pairwiseFunctor, traversalSelectorInfo, datalayout, newton3);

    checkTraversal(traversalOption, std::move(traversal));
  }

  MTriwiseFunctor triwiseFunctor;
  for (const auto &traversalOption : autopas::TraversalOption::getAllTriwiseOptions()) {
    auto traversal = autopas::TraversalSelector::generateTraversal<FPCell, MTriwiseFunctor>(
        traversalOption, triwiseFunctor, traversalSelectorInfo, autopas::DataLayoutOption::aos, false);

    checkTraversal(traversalOption, std::move(traversal));
  }
}
