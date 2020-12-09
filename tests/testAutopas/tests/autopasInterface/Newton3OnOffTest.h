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
#include "autopas/options/TraversalOption.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

namespace Newton3OnOffTest {

/**
 * Test to check if newton3 and non-newton3 work as expected
 */
class Newton3OnOffTest
    : public AutoPasTestBase,
      public ::testing::WithParamInterface<
          std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption>> {
 public:
  Newton3OnOffTest() : mockFunctor() {}

  void SetUp() override {}

  void TearDown() override {}

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {10.0, 10.0, 10.0}; }

  static double getCutoff() { return 1.0; }
  static double getCellSizeFactor() { return 1.0; }
  static double getVerletSkin() { return 0.0; }
  static int getClusterSize() { return 4; }

  void countFunctorCalls(autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
                         autopas::DataLayoutOption dataLayout);

  template <class Container, class Traversal>
  void iterate(Container container, Traversal traversal);

  MockFunctor<Particle> mockFunctor;

  /**
   * Determines how often the functor is called for single cells and pairs of cells und run additional checks.
   * @tparam useNewton3 Enables or disables newton3.
   * @tparam Container Type of container.
   * @tparam Traversal Type of traversal.
   * @param dataLayout Data layout.
   * @param container Container.
   * @param traversalOption Traversal option.
   * @return [#calls single cell, #calls pair of cells]
   */
  template <bool useNewton3, class Container, class Traversal>
  std::pair<size_t, size_t> eval(autopas::DataLayoutOption dataLayout, Container &container, Traversal traversalOption);

  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto inputTuple = static_cast<ParamType>(info.param);

      auto [containerOption, traversalOption, dataLayoutOption] = inputTuple;

      auto retStr =
          containerOption.to_string() + "_" + traversalOption.to_string() + "_" + dataLayoutOption.to_string();
      // replace all '-' with '_', otherwise the test name is invalid
      std::replace(retStr.begin(), retStr.end(), '-', '_');
      return retStr;
    }
  };
};

}  // end namespace Newton3OnOffTest
