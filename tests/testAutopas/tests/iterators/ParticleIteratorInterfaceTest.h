/**
 * @file ParticleIteratorInterfaceTest.h
 * @author seckler
 * @date 22.07.19
 */

#pragma once

#include <gtest/gtest.h>

#include <tuple>

#include "autopas/AutoPas.h"
#include "testingHelpers/commonTypedefs.h"

using testingTuple =
    std::tuple<autopas::ContainerOption, double /*cell size factor*/, bool /*regionIterator (true) or regular (false)*/,
               bool /*testConstIterators*/, bool /*priorForceCalc*/, autopas::IteratorBehavior>;

class ParticleIteratorInterfaceTest : public testing::Test, public ::testing::WithParamInterface<testingTuple> {
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
   * Initialize the given AutoPas object with the default
   * @tparam AutoPasT
   * @param autoPas
   * @return tuple {haloBoxMin, haloBoxMax}
   */
  template <typename AutoPasT>
  auto defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption, double cellSizeFactor);

  /**
   * Apply an iterator, track what particle IDs are found and compare this to a vector of expected IDs
   * @tparam AutoPasT
   * @tparam IteratorT
   * @param autopas
   * @param iterator
   * @param particleIDsExpected
   */
  template <class AutoPasT, class FgetIter>
  void findParticles(AutoPasT &autopas, FgetIter getIter, const std::vector<size_t> &particleIDsExpected);

  /**
   * Inserts particles around all corners of the given AutoPas object at critical distances.
   * @tparam AutoPasT
   * @param autoPas
   * @return Tuple of two vectors containing IDs of added particles. First for owned, second for halo particles.
   */
  template <class AutoPasT>
  auto fillContainerAroundBoundary(AutoPasT &autoPas);

  /**
   * Creats a grid of particles in the given AutoPas object.
   * Grid width is 3 *( boxLength / ((cutoff + skin) * cellSizeFactor) ).
   * The lower corner of the grid is offset from boxMin by half the grid width in every dimension.
   * This way there should be one particle in every third Linked Cells cell.
   * @tparam AutoPasT
   * @param autoPas
   * @return Vector of all particle IDs added.
   */
  template <class AutoPasT>
  auto fillContainerWithGrid(AutoPasT &autoPas);

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
  auto deleteParticles(AutoPasT &autopas, F predicate, bool useRegionIterator, autopas::IteratorBehavior behavior);

  /**
   * Creates a function to instantiate an iterator with the given properties and passes this function on to fun.
   * This is necessary so that fun can decide for itself if it wants Iterators to be created in an OpenMP region or not.
   * @tparam AutoPasT
   * @tparam F f(AutoPas, Iterator)
   * @param useRegionIterator
   * @param useConstIterator
   * @param behavior
   * @param autoPas
   * @param fun Function taking the AutoPas object and the generated iterator.
   */
  template <class AutoPasT, class F>
  void provideIterator(bool useRegionIterator, bool useConstIterator, autopas::IteratorBehavior behavior,
                       AutoPasT &autoPas, F fun);

  /**
   * Same as provideIterator() but with const const indicator.
   * @tparam useConstIterator
   * @tparam AutoPasT
   * @tparam F f(AutoPas, Iterator)
   * @param useRegionIterator
   * @param behavior
   * @param autoPas
   * @param fun Function taking the AutoPas object and the generated iterator.
   */
  template <bool useConstIterator, class AutoPasT, class F>
  void provideIterator(bool useRegionIterator, autopas::IteratorBehavior behavior, AutoPasT &autoPas, F fun);

  template <bool testConstIterators>
  void testAdditionAndIteration(autopas::ContainerOption containerOption, double cellSizeOption, bool priorForceCalc);

  template <bool testConstIterators>
  void testOpenMPIterators(autopas::ContainerOption containerOption, double cellSizeFactor,
                           autopas::IteratorBehavior behavior, bool testRegionIterators, bool priorForceCalc);

  /**
   * The main idea of this test is to check whether deletion through a RegionIterator doesn't do stupid things.
   * Step by step:
   * - Fill AutoPas object with 200 random owned particles (might be in halo)
   * - Do or don't do a force calculation
   * - Save all generated particles in a vector
   * - Use region iterator to delete all particles from a region and save what was deleted
   * - Use regular Iterator to check everything that should be there is there
   */
  template <bool testConstIterators>
  void testRegionIteratorDeletion(autopas::ContainerOption containerOption, double cellSizeFactor, bool priorForceCalc);

  /**
   * Tests the equivalence of the range-based for loop with the normal for loop using isValid
   * @param containerOption
   */
  template <bool testConstIterators>
  void testRangeBasedIterator(autopas::ContainerOption containerOption, double cellSizeOption, bool priorForceCalc);
};