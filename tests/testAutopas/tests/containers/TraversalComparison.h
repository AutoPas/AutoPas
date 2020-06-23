/**
 * @file TraversalComparison.h
 * @author humig
 * @date 12.07.19
 */

#pragma once

#include <gtest/gtest.h>

#include <cstdlib>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/commonTypedefs.h"

enum DeletionPosition {
  // We have chosen the values explicitly, s.t., this enum can be used using bit manipulation, i.e., beforeAndAfterLists
  // enables both bits for beforeLists and afterLists.
  never = 0,
  beforeLists = 1,
  afterLists = 2,
  beforeAndAfterLists = 3
};

using TestingTuple = std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption,
                                autopas::Newton3Option, size_t /*numParticles*/, size_t /*numHaloParticles*/,
                                std::array<double, 3> /*boxMaxVec*/, double /*cellSizeFactor*/, bool /*doSlightShift*/,
                                DeletionPosition /*particleDeletionPosition*/, bool /*globals*/>;
/**
 * The tests in this class compare the calculated forces from all aos and soa traversals with a reference result.
 */
class TraversalComparison : public AutoPasTestBase, public ::testing::WithParamInterface<TestingTuple> {
 public:
  using mykey_t = std::tuple<size_t,                 // numParticles
                             size_t,                 // numHaloParticles
                             std::array<double, 3>,  // boxMax
                             bool,                   // doSlightShift
                             DeletionPosition,       // particleDeletionPosition
                             bool                    // globals
                             >;

  /**
   * Struct to hold global values
   */
  struct Globals {
    double upot{};
    double virial{};
  };

  static void generateReference(mykey_t key);

  static auto getTestParams();

 protected:
  template <class ContainerPtrType>
  static void executeShift(ContainerPtrType containerPtr, double magnitude, size_t totalNumParticles);

  template <bool globals>
  static std::tuple<std::vector<std::array<double, 3>>, Globals> calculateForces(
      autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
      autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option, size_t numMolecules,
      size_t numHaloMolecules, std::array<double, 3> boxMax, double cellSizeFactor, bool doSlightShift,
      DeletionPosition particleDeletionPosition);

  static constexpr std::array<double, 3> _boxMin{0, 0, 0};
  static constexpr double _cutoff{1.};

  static constexpr double _eps{1.};
  static constexpr double _sig{1.};

  static inline std::map<mykey_t, std::vector<std::array<double, 3>>> _forcesReference{};
  static inline std::map<mykey_t, Globals> _globalValuesReference{};
};
