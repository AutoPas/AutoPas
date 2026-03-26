/**
 * @file Newton3OnOffTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/TraversalOption.h"
#include "mocks/MockPairwiseFunctor.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Test to check if newton3 and non-newton3 work as expected
 */
class Newton3OnOffTest : public AutoPasTestBase, public ::testing::WithParamInterface<autopas::Configuration> {
 public:
  Newton3OnOffTest() {}

  void SetUp() override {}

  static std::array<double, 3> getBoxMin() { return {0.0, 0.0, 0.0}; }
  static std::array<double, 3> getBoxMax() { return {10.0, 10.0, 10.0}; }

  static double getCutoff() { return 1.0; }
  static double getVerletSkin() { return 0.0; }
  static int getClusterSize() { return 4; }
  static int getSortingThreshold() { return 8; }

  /**
   * Count the number of Functor calls with and without newton 3 and compare.
   *
   * For AoS,
   * - for pairwise, there should be twice as many functor calls without newton3 as with
   * - for triwise, there should be thrice as many functor calls without newton3 as with
   *
   * For SoA, regardless of interaction type,
   * - for SoASingle, there should be the same number of functor calls with and without newton3
   * - for SoAPair, there should be twice as many functor calls without newton3 as with
   * - for SoATriplet, there should be thrice as many functor calls without newton3 as with
   *
   * @tparam Functor_T (Mock) functor type
   * @param config
   */
  template <typename Functor_T>
  void countFunctorCalls(autopas::Configuration config);

  /**
   * Mimics a simulation iteration (rebuilds neighbor lists, the computes interactions)
   * @tparam Container_T container type
   * @tparam Traversal_T traversal type
   * @param container the actual container.
   * @param traversal the actual traversal.
   */
  template <class Container_T, class Traversal_T>
  void iterate(Container_T &container, Traversal_T traversal);

  /**
   * Determines how often the functor is called for single cells and pairs of cells und run additional checks.
   * @tparam Functor_T Type of functor.
   * @tparam Container_T Type of container.
   * @param config Configuration. Should include newton3 disabled, regardless of useNewton3.
   * @param useNewton3 If true, a newton3 enabled configuration is created.
   * @param container The actual container.
   * @param mockFunctor Functor that is used for evaluation.
   * @return [#calls single cell, #calls pair of cells, #calls triplet of cells (in case of 3b)]
   */
  template <class Functor_T, class Container_T>
  std::tuple<size_t, size_t, size_t> eval(autopas::Configuration config, bool useNewton3, Container_T &container,
                                          Functor_T &mockFunctor);

  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto config = info.param;

      auto retStr = config.toShortString(false, true);
      return retStr;
    }
  };
};
