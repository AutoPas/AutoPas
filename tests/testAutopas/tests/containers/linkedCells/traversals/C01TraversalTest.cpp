/**
 * @file C01TraversalTest.cpp
 * @author S. Seckler
 * @date 10.01.2019
 */

#include "C01TraversalTest.h"

#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "testingHelpers/commonTypedefs.h"

namespace C01TraversalTest {

// Place to implement special test cases, which only apply to C01 Traversal

TEST_F(C01TraversalTest, testIsApplicable) {
  const std::array<unsigned long, 3> dims({1ul, 1ul, 1ul});
  MFunctor functor;

  autopas::LCC01Traversal<FPCell, MFunctor, autopas::DataLayoutOption::soa, false, false> c01T_N3off(dims, &functor, 1.,
                                                                                                     {1., 1., 1.});
  EXPECT_EQ(c01T_N3off.isApplicable(), true);

  autopas::LCC01Traversal<FPCell, MFunctor, autopas::DataLayoutOption::soa, true, false> c01T_N3on(dims, &functor, 1.,
                                                                                                   {1., 1., 1.});
  EXPECT_EQ(c01T_N3on.isApplicable(), false);

  autopas::LCC01Traversal<FPCell, MFunctor, autopas::DataLayoutOption::soa, false, true> c01T_N3off_combineSoA(
      dims, &functor, 1., {1., 1., 1.});
  EXPECT_EQ(c01T_N3off_combineSoA.isApplicable(), true);

  autopas::LCC01Traversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, false, true> c01T_N3off_combineSoA_AoS(
      dims, &functor, 1., {1., 1., 1.});
  EXPECT_EQ(c01T_N3off_combineSoA_AoS.isApplicable(), false);

  autopas::LCC01Traversal<FPCell, MFunctor, autopas::DataLayoutOption::soa, true, true> c01T_N3on_combineSoA(
      dims, &functor, 1., {1., 1., 1.});
  EXPECT_EQ(c01T_N3on_combineSoA.isApplicable(), false);
}
}  // end namespace C01TraversalTest
