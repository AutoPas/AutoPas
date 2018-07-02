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

  autopas::TraversalSelector<FPCell> traversalSelectorC08({30, 30, 30}, optionVectorC08);
  autopas::TraversalSelector<FPCell> traversalSelectorSlice({30, 30, 30}, optionVectorSliced);

  auto traversalC08 = traversalSelectorC08.getOptimalTraversal(cf);
  auto traversalSlice = traversalSelectorSlice.getOptimalTraversal(cf);

  // check that traversals are of the expected type
  ASSERT_TRUE((dynamic_cast<autopas::C08Traversal<FPCell, CellFunctorAoSN3> *>(traversalC08.get())));
  ASSERT_TRUE((dynamic_cast<autopas::SlicedTraversal<FPCell, CellFunctorAoSN3> *>(traversalSlice.get())));
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