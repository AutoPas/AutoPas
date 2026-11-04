/**
 * @file TraversalSelectorTest.cpp
 * @author F. Gratl
 * @date 21.06.18
 */

#include "TraversalSelectorTest.h"

#include <type_traits>

#include "autopas/containers/cellTraversals/C01BasedTraversal.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCAllCellsNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/neighborLists/VLCCellPairNeighborList.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCC01Traversal.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::Return;

using C01TraversalBase = autopas::C01BasedTraversal<FPCell, MPairwiseFunctor>;
static_assert(
    std::is_base_of_v<C01TraversalBase, autopas::VLCC01Traversal<FPCell, MPairwiseFunctor,
                                                                 autopas::VLCAllCellsNeighborList<ParticleFP64>>>,
    "vlc_c01 must use the default collapse depth.");
static_assert(
    std::is_base_of_v<C01TraversalBase, autopas::VLCC01Traversal<FPCell, MPairwiseFunctor,
                                                                 autopas::VLCCellPairNeighborList<ParticleFP64>>>,
    "vlp_c01 must use the default collapse depth.");

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
    auto onlySoATraversals = autopas::compatibleTraversals::allTraversalsSupportingOnlySoA();
    if (std::find(onlySoATraversals.begin(), onlySoATraversals.end(), traversalOption) != onlySoATraversals.end()) {
      datalayout = autopas::DataLayoutOption::soa;  // these traversals only work with SoA
    }
    auto onlyN3Traversals = autopas::compatibleTraversals::allTraversalsSupportingOnlyNewton3Enabled();
    if (std::find(onlyN3Traversals.begin(), onlyN3Traversals.end(), traversalOption) != onlyN3Traversals.end()) {
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
