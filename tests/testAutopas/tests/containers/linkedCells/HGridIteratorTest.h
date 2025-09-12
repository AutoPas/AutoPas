/**
 * @file HGridIteratorTest.h
 * @author atacann
 * @date 15.01.2024
 */

#pragma once

#include <gtest/gtest.h>

#include <tuple>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/IteratorBehavior.h"

class HGridIteratorTestBase : public testing::Test {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto [containerOption, cellSizeFactor, useRegionIterator, testConstIterators, priorForceCalc, behavior] =
          static_cast<ParamType>(info.param);
      std::string str;
      str += containerOption.to_string() + "_";
      str += std::string{useRegionIterator ? "Region" : ""} + "Iterator_";
      str += std::string{"cellSizeFactor"} + std::to_string(cellSizeFactor);
      str += testConstIterators ? "_const" : "_nonConst";
      str += priorForceCalc ? "_priorForceCalc" : "_noPriorForceCalc";
      str += "_" + behavior.to_string();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };

  /**
   * Initialize the given AutoPas object with the default values for this test class.
   * @tparam AutoPasT
   * @param autoPas
   * @return tuple {haloBoxMin, haloBoxMax}
   */
  template <typename AutoPasT>
  auto defaultInit(AutoPasT &autoPas, const autopas::ContainerOption &containerOption, double cellSizeFactor);

  /**
   * Deletes all particles whose ID matches the given Predicate.
   * @tparam constIter
   * @tparam AutoPasT
   * @tparam F
   * @param autopas
   * @param predicate
   * @param useRegionIterator
   * @param behavior
   * @return
   */
  template <bool constIter, class AutoPasT, class F>
  auto deleteParticles(AutoPasT &autopas, F predicate, bool useRegionIterator,
                       const autopas::IteratorBehavior &behavior);
};

using testingTuple =
    std::tuple<autopas::ContainerOption, double /*cell size factor*/, bool /*regionIterator (true) or regular (false)*/,
               bool /*testConstIterators*/, bool /*priorForceCalc*/, autopas::IteratorBehavior>;
class HGridIteratorTest : public HGridIteratorTestBase, public ::testing::WithParamInterface<testingTuple> {};

class HGridIteratorTestNonConst : public HGridIteratorTestBase, public ::testing::WithParamInterface<testingTuple> {};

class HGridIteratorTestNonConstOwned : public HGridIteratorTestBase,
                                       public ::testing::WithParamInterface<testingTuple> {};
