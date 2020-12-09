/**
 * @file TraversalSelectorTest.cpp
 * @author F. Gratl
 * @date 21.06.18
 */

#include "TraversalSelectorTest.h"

namespace TraversalSelectorTest {

using ::testing::Return;

/**
 * Check if the only allowed traversal is returned
 */
TEST_F(TraversalSelectorTest, testSelectAndGetCurrentTraversal) {
  MFunctor functor;

  // this should be high enough so that sliced is still valid for the current processors thread count.
  constexpr size_t domainSize = 900;
  autopas::TraversalSelectorInfo traversalSelectorInfo({domainSize, domainSize, domainSize}, 1., {1., 1., 1.}, 32);

  for (const auto &traversalOption : autopas::TraversalOption::getAllOptions()) {
    auto traversal =
        autopas::TraversalSelector<FPCell>::template generateTraversal<MFunctor, autopas::DataLayoutOption::aos, false>(
            traversalOption, functor, traversalSelectorInfo);

    // check that traversals are of the expected type
    EXPECT_EQ(traversalOption, traversal->getTraversalType())
        << "Is the domain size large enough for the processors' thread count?";
  }
}

}  // end namespace TraversalSelectorTest
