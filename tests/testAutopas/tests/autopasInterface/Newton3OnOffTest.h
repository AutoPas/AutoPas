/**
 * @file Newton3OnOffTest.h
 * @author seckler
 * @date 18.04.18
 */

#pragma once

#include <gtest/gtest.h>
#include <testingHelpers/commonTypedefs.h>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "autopas/sph/autopassph.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/RandomGenerator.h"

/**
 * Test to check if newton3 and non-newton3 work as expected
 */
class Newton3OnOffTest : public AutoPasTestBase,
                         public ::testing::WithParamInterface<std::tuple<std::string, std::string>> {
 public:
  Newton3OnOffTest() : mockFunctor() {}

  void SetUp() override {}

  void TearDown() override {}

  std::array<double, 3> getBoxMin() const { return {0.0, 0.0, 0.0}; }

  std::array<double, 3> getBoxMax() const { return {3.0, 3.0, 3.0}; }

  double getCutoff() const { return 1.0; }
  double getCellSizeFactor() const { return 1.0; }
  double getVerletSkin() const { return 0.0; }

  void countFunctorCalls(autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
                         autopas::DataLayoutOption dataLayout);

  template <class ParticleFunctor, class Container, class Traversal>
  void iterate(Container container, Traversal traversal, autopas::DataLayoutOption dataLayout,
               autopas::Newton3Option newton3, ParticleFunctor *f);

  MockFunctor<Particle, FPCell> mockFunctor;

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
};
