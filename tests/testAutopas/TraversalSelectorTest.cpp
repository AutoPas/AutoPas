/**
 * @file TraversalSelectorTest.cpp
 * @author F. Gratl
 * @date 21.06.18
 */
#include "TraversalSelectorTest.h"

/**
 * Check if the only allowed traversal is returned
 */
TEST(TraversalSelectorTest, testGetOptimalTraversalOneOption) {
  MockFunctor<autopas::Particle, FPCell> functor;
  CellFunctorAoSN3 cf(&functor);

  std::vector<autopas::TraversalOptions> optionVectorC08 = {autopas::TraversalOptions::c08};
  std::vector<autopas::TraversalOptions> optionVectorSliced = {autopas::TraversalOptions::sliced};

  // this should be high enough so that sliced is still valid for the current processors thread count.
  constexpr size_t domainSize = 1000;

  autopas::TraversalSelector<FPCell> traversalSelectorC08({domainSize, domainSize, domainSize}, optionVectorC08);
  autopas::TraversalSelector<FPCell> traversalSelectorSlice({domainSize, domainSize, domainSize}, optionVectorSliced);

  auto traversalC08 = traversalSelectorC08.getOptimalTraversal(cf);
  auto traversalSlice = traversalSelectorSlice.getOptimalTraversal(cf);

  // check that traversals are of the expected type
  ASSERT_TRUE((dynamic_cast<autopas::C08Traversal<FPCell, CellFunctorAoSN3> *>(traversalC08.get())));
  ASSERT_TRUE((dynamic_cast<autopas::SlicedTraversal<FPCell, CellFunctorAoSN3> *>(traversalSlice.get())))
      << "Is the domain size large enough for the processors' thread count?";
}

TEST(TraversalSelectorTest, testGetOptimalTraversalBadFirstOption) {
  MockFunctor<autopas::Particle, FPCell> functor;
  CellFunctorAoSN3 cf(&functor);

  std::vector<autopas::TraversalOptions> optionVector = {autopas::TraversalOptions::sliced,
                                                         autopas::TraversalOptions::c08};

  autopas::TraversalSelector<FPCell> traversalSelectorC08({1, 1, 1}, optionVector);

  auto traversal = traversalSelectorC08.getOptimalTraversal(cf);

  // check that traversals are of the expected type
  ASSERT_TRUE((dynamic_cast<autopas::C08Traversal<FPCell, CellFunctorAoSN3> *>(traversal.get())));
}