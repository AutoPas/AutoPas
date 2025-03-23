/**
 * @file HGridTraversalComparison.h
 * @author atacann
 * @date 15.01.2024
 */

// Mostly copied from TraversalComparison.h, minor changes to use different molecule types and test HGrid only

#pragma once

#include <gtest/gtest.h>

#include <cstdlib>

#include "AutoPasTestBase.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/AxilrodTellerFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

enum DeletionPosition {
  // We have chosen the values explicitly, s.t., this enum can be used using bit manipulation, i.e., beforeAndAfterLists
  // enables both bits for beforeLists and afterLists.
  never = 0,
  beforeLists = 1,
  afterLists = 2,
  beforeAndAfterLists = 3
};

using TestingTuple =
    std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption, autopas::Newton3Option,
               size_t /*numParticles*/, size_t /*numHaloParticles*/, std::array<double, 3> /*boxMaxVec*/,
               double /*cellSizeFactor*/, bool /*doSlightShift*/, DeletionPosition /*particleDeletionPosition*/,
               bool /*globals*/, autopas::InteractionTypeOption>;

/**
 * The tests in this class compare the calculated forces from all aos and soa traversals with a reference result.
 */
class HGridTraversalComparison : public AutoPasTestBase, public ::testing::WithParamInterface<TestingTuple> {
 public:
  using mykey_t = std::tuple<size_t,                                  // numParticles
                             size_t,                                  // numHaloParticles
                             std::array<double, 3>,                   // boxMax
                             bool,                                    // doSlightShift
                             DeletionPosition,                        // particleDeletionPosition
                             bool,                                    // globals
                             autopas::InteractionTypeOption::Value>;  // interaction type

  /**
   * Struct to hold global values
   */
  struct Globals {
    double upot{};
    double virial{};
  };

  /**
   * Struct to hold parameters that might differ for different interaction types
   */
  struct TraversalTestParams {
    double deletionPercentage;
    std::vector<unsigned long> numParticles;
    std::vector<unsigned long> numHaloParticles;
    std::vector<std::array<double, 3>> boxMax;
    std::vector<double> cellSizeFactors;
  };

  template <bool globals>
  static void generateReference(mykey_t key);

  static auto getTestParams();

  static std::array<double, 3> randomShift(double magnitude, std::mt19937 &generator);

  static std::unordered_map<autopas::InteractionTypeOption::Value, TraversalTestParams> params;

 protected:
  template <class ContainerType>
  static void executeShift(ContainerType &containerPtr, double magnitude, size_t numTotalParticles);

  template <typename ContainerT>
  static void markSomeParticlesAsDeleted(ContainerT &container, size_t numTotalParticles, unsigned seed,
                                         autopas::InteractionTypeOption interactionT);

  template <bool globals>
  static std::tuple<std::vector<std::array<double, 3>>, Globals> calculateForces(
      autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
      autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, double cellSizeFactor,
      mykey_t key, bool useSorting);

  template <typename Functor, bool globals>
  static std::tuple<std::vector<std::array<double, 3>>, Globals> calculateForcesImpl(
      Functor functor, autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
      autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, double cellSizeFactor,
      mykey_t key, bool useSorting);

  static constexpr std::array<double, 3> _boxMin{0, 0, 0};
  static constexpr double _cutoff{1.};

  static std::shared_ptr<ParticlePropertiesLibrary<>> _particlePropertiesLibrary;

  static inline std::map<mykey_t, std::vector<std::array<double, 3>>> _forcesReference{};
  static inline std::map<mykey_t, Globals> _globalValuesReference{};
};
