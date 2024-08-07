/**
 * @file C01TraversalTest.cpp
 * @author S. Seckler
 * @date 10.01.2019
 */

#include "C01TraversalTest.h"

#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "autopas/utils/checkFunctorType.h"
#include "testingHelpers/commonTypedefs.h"

// Place to implement special test cases, which only apply to C01 Traversal

TEST_F(C01TraversalTest, testIsApplicable) {
  const std::array<unsigned long, 3> dims({1ul, 1ul, 1ul});

  /// Pairwise traversal tests
  MPairwiseFunctor pairwiseFunctor;

  autopas::LCC01Traversal<FPCell, MPairwiseFunctor, false> c01T_N3off(dims, &pairwiseFunctor, 1., {1., 1., 1.},
                                                                      autopas::DataLayoutOption::soa, false);
  EXPECT_EQ(c01T_N3off.isApplicable(), true);

  autopas::LCC01Traversal<FPCell, MPairwiseFunctor, false> c01T_N3on(dims, &pairwiseFunctor, 1., {1., 1., 1.},
                                                                     autopas::DataLayoutOption::soa, true);
  EXPECT_EQ(c01T_N3on.isApplicable(), false);

  autopas::LCC01Traversal<FPCell, MPairwiseFunctor, true> c01T_N3off_combineSoA(
      dims, &pairwiseFunctor, 1., {1., 1., 1.}, autopas::DataLayoutOption::soa, false);
  EXPECT_EQ(c01T_N3off_combineSoA.isApplicable(), true);

  autopas::LCC01Traversal<FPCell, MPairwiseFunctor, true> c01T_N3off_combineSoA_AoS(
      dims, &pairwiseFunctor, 1., {1., 1., 1.}, autopas::DataLayoutOption::aos, false);
  EXPECT_EQ(c01T_N3off_combineSoA_AoS.isApplicable(), false);

  autopas::LCC01Traversal<FPCell, MPairwiseFunctor, true> c01T_N3on_combineSoA(dims, &pairwiseFunctor, 1., {1., 1., 1.},
                                                                               autopas::DataLayoutOption::soa, true);
  EXPECT_EQ(c01T_N3on_combineSoA.isApplicable(), false);

  /// Triwise traversal tests
  MTriwiseFunctor triwiseFunctor;

  autopas::LCC01Traversal<FPCell, MTriwiseFunctor> c01T_N3off_3B(dims, &triwiseFunctor, 1., {1., 1., 1.},
                                                                 autopas::DataLayoutOption::aos, false);
  EXPECT_EQ(c01T_N3off_3B.isApplicable(), true);

  autopas::LCC01Traversal<FPCell, MTriwiseFunctor> c01T_N3on_3B(dims, &triwiseFunctor, 1., {1., 1., 1.},
                                                                autopas::DataLayoutOption::aos, true);
  EXPECT_EQ(c01T_N3on_3B.isApplicable(), false);
}