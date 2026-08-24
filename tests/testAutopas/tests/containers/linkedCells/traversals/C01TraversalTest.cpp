/**
 * @file C01TraversalTest.cpp
 * @author S. Seckler
 * @date 10.01.2019
 */

#include "C01TraversalTest.h"

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LCC01Traversal.h"
#include "autopas/tuning/Configuration.h"
#include "testingHelpers/commonTypedefs.h"

// Place to implement special test cases, which only apply to C01 Traversal

TEST_F(C01TraversalTest, testHasCompatibleValues) {
  using namespace autopas;

  /// Pairwise traversal tests
  Configuration c01T_N3off(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01, LoadEstimatorOption::none,
                           DataLayoutOption::soa, Newton3Option::disabled, InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3off.hasCompatibleValues(), true);

  Configuration c01T_N3off_AoS(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01,
                               LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::disabled,
                               InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3off_AoS.hasCompatibleValues(), true);

  Configuration c01T_N3on(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01, LoadEstimatorOption::none,
                          DataLayoutOption::soa, Newton3Option::enabled, InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3on.hasCompatibleValues(), false);

  Configuration c01T_N3on_AoS(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01,
                              LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::enabled,
                              InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3on_AoS.hasCompatibleValues(), false);

  // Combined SoA Tests

  Configuration c01T_N3off_combineSoA(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA,
                                      LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::disabled,
                                      InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3off_combineSoA.hasCompatibleValues(), true);

  Configuration c01T_N3off_combineSoA_AoS(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA,
                                          LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::disabled,
                                          InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3off_combineSoA_AoS.hasCompatibleValues(), false);

  Configuration c01T_N3on_combineSoA(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA,
                                     LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::enabled,
                                     InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3on_combineSoA.hasCompatibleValues(), false);

  Configuration c01T_N3on_combineSoA_AoS(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA,
                                         LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::enabled,
                                         InteractionTypeOption::pairwise);
  EXPECT_EQ(c01T_N3on_combineSoA_AoS.hasCompatibleValues(), false);

  /// Triwise traversal tests
  Configuration c01T_N3off_3B(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01,
                              LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::disabled,
                              InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3off_3B.hasCompatibleValues(), true);

  Configuration c01T_N3off_AoS_3B(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01,
                                  LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::disabled,
                                  InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3off_AoS_3B.hasCompatibleValues(), true);

  Configuration c01T_N3on_3B(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01,
                             LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::enabled,
                             InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3on_3B.hasCompatibleValues(), false);

  Configuration c01T_N3on_AoS_3B(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01,
                                 LoadEstimatorOption::none, DataLayoutOption::aos, Newton3Option::enabled,
                                 InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3on_AoS_3B.hasCompatibleValues(), false);

  // Combined SoA Tests

  Configuration c01T_N3off_combineSoA_3B(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA,
                                         LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::disabled,
                                         InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3off_combineSoA_3B.hasCompatibleValues(), false);

  Configuration c01T_N3off_combineSoA_AoS_3B(
      ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA, LoadEstimatorOption::none,
      DataLayoutOption::aos, Newton3Option::disabled, InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3off_combineSoA_AoS_3B.hasCompatibleValues(), false);

  Configuration c01T_N3on_combineSoA_3B(ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA,
                                        LoadEstimatorOption::none, DataLayoutOption::soa, Newton3Option::enabled,
                                        InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3on_combineSoA_3B.hasCompatibleValues(), false);

  Configuration c01T_N3on_combineSoA_AoS_3B(
      ContainerOption::Value::linkedCells, 1., TraversalOption::lc_c01_combined_SoA, LoadEstimatorOption::none,
      DataLayoutOption::aos, Newton3Option::enabled, InteractionTypeOption::triwise);
  EXPECT_EQ(c01T_N3on_combineSoA_AoS_3B.hasCompatibleValues(), false);
}