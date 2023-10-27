/**
 * @file C01TraversalTest3B.cpp
 * @author muehlhaeusser
 * @date 26.10.2023
 */

#include "C01TraversalTest3B.h"

#include "autopas/containers/linkedCells/traversals/LCC01Traversal3B.h"
#include "testingHelpers/commonTypedefs.h"

// Place to implement special test cases, which only apply to C01 Traversal

TEST_F(C01TraversalTest3B, testIsApplicable) {
  const std::array<unsigned long, 3> dims({1ul, 1ul, 1ul});
  MTriwiseFunctor functor;

  autopas::LCC01Traversal3B<FPCell, MTriwiseFunctor, autopas::DataLayoutOption::aos, false> c01T_N3off(
      dims, &functor, 1., {1., 1., 1.});
  EXPECT_EQ(c01T_N3off.isApplicable(), true);

  autopas::LCC01Traversal3B<FPCell, MTriwiseFunctor, autopas::DataLayoutOption::aos, true> c01T_N3on(dims, &functor, 1.,
                                                                                                     {1., 1., 1.});
  EXPECT_EQ(c01T_N3on.isApplicable(), false);
}