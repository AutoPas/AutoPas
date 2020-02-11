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
using TestingTuple = std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption,
                                autopas::Newton3Option, size_t /*numParticles*/, std::array<double, 3> /*boxMaxVec*/>;
/**
 * The tests in this class compare the calculated forces from all aos and soa traversals with a reference result.
 */
class TraversalComparison : public AutoPasTestBase, public ::testing::WithParamInterface<TestingTuple> {
 public:
  static void SetUpTestSuite();

  static auto getTestParams();

 protected:
  static std::vector<std::array<double, 3>> calculateForces(autopas::ContainerOption containerOption,
                                                            autopas::TraversalOption traversalOption,
                                                            autopas::DataLayoutOption dataLayoutOption,
                                                            autopas::Newton3Option newton3Option,
                                                            unsigned long numMolecules, std::array<double, 3> boxMax);

  static inline std::array<double, 3> _boxMin{0, 0, 0};
  static inline std::array<std::array<double, 3>, 2> _boxMaxVector{{{3, 3, 3}, {10, 10, 10}}};
  static inline double _cutoff{1.};

  static inline double _eps{1.};
  static inline double _sig{1.};

  static inline std::map<std::pair<int, std::array<double, 3>>, std::vector<std::array<double, 3>>> _forcesReference{};

  static inline std::array<int, 3> _numParticlesVector{100, 1000, 2000};
};
