/**
 * @file NeighborIdentificationFunctorTest.h
 * @author S. Newcome
 * @date 16.07.2026
 */

#pragma once

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "AutoPasTestBase.h"
#include "autopas/particles/OwnershipState.h"

/**
 * Base fixture for all NeighborIdentificationFunctor tests.
 */
class NeighborIdentificationFunctorTest : public AutoPasTestBase {};

/**
 * Per-particle neighbor lists given as particle IDs, indexed by particle ID.
 */
using ExpectedLists = std::vector<std::vector<size_t>>;

/**
 * Name generator for parameterized test cases. Uses the name member of the parameter struct as the test's name suffix.
 * Parameters combined via ::testing::Combine are handled by the tuple overload, which concatenates the names of the
 * combined components.
 */
struct ParamNameGenerator {
  /**
   * Returns the name of the given test case.
   * @tparam ParamType Type of the parameter struct. Must have a name member.
   * @param info Test parameter info provided by gtest.
   * @return The parameter's name.
   */
  template <class ParamType>
  std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
    return info.param.name;
  }

  /**
   * Returns the name of the given test case, where the test is a combination of parameter components.
   * @tparam ParamTypes Types of the parameter structs. Each must have a name member.
   * @param info Test parameter info provided by gtest.
   * @return The parameters' concatenated names.
   */
  template <class... ParamTypes>
  std::string operator()(const testing::TestParamInfo<std::tuple<ParamTypes...>> &info) const {
    return std::apply([](const auto &...params) { return (params.name + ...); }, info.param);
  }
};

/**
 * Describes the newton3 usage of one test case - whether newton3 is used in list generation and whether newton3 lists
 * are themselves generated.
 */
struct N3Mode {
  /**
   * Name of the mode.
   */
  std::string name;

  /**
   * Whether newton3 lists are gathered.. In the full-flow tests this also determines whether newton3 is used when
   * calculating the forces using the generated lists.
   */
  bool gatherN3Lists;

  /**
   * Whether newton3 is used for list generation by InteractionListGeneratorFunctor.
   */
  bool n3Used;
};

/**
 * Describes the ownership of the particles of one test case.
 */
struct OwnershipVariant {
  /**
   * Name of the variant, used as part of the test case's name.
   */
  std::string name;

  /**
   * OwnershipStates assigned to the particles in a round-robin. A single state gives every particle that state, two
   * states alternate between them, etc...
   */
  std::vector<autopas::OwnershipState> states;
};

/**
 * Parameters describing an interaction list generation test case. Use for all functors except SoAFunctorPair (see
 * PairListGenerationParams)
 */
struct ListGenerationParams {
  /**
   * Name of the test case.
   */
  std::string name;

  /**
   * Position of each particle. One particle is created per entry.
   */
  std::vector<std::array<double, 3>> positions;

  /**
   * OwnershipState of each particle. Must have the same size as positions.
   */
  std::vector<autopas::OwnershipState> ownerships;

  /**
   * The newton3 usage of this test case.
   */
  N3Mode n3Mode;

  /**
   * Expected neighbor list of each particle, given as ascending particle IDs. One entry per particle.
   */
  ExpectedLists expectedNeighbors;

  /**
   * Full Verlet list of each particle, assuming N3 not used (pruned as needed). Only used by the SoAFunctorVerlet
   * tests.
   */
  ExpectedLists verletLists;
};

/**
 * Describes one SoAFunctorPair test case operating on two cells.
 */
struct PairListGenerationParams {
  /**
   * Name of the test case.
   */
  std::string name;

  /**
   * Position of each particle of the first cell.
   */
  std::vector<std::array<double, 3>> positions1;

  /**
   * Position of each particle of the second cell.
   */
  std::vector<std::array<double, 3>> positions2;

  /**
   * OwnershipState of each particle of the first cell. Must have the same size as positions1.
   */
  std::vector<autopas::OwnershipState> ownerships1;

  /**
   * OwnershipState of each particle of the second cell. Must have the same size as positions2.
   */
  std::vector<autopas::OwnershipState> ownerships2;

  /**
   * The newton3 usage of this test case.
   */
  N3Mode n3Mode;

  /**
   * Expected neighbor list of each particle of the first cell, given as ascending particle IDs.
   */
  ExpectedLists expectedNeighbors1;

  /**
   * Expected neighbor list of each particle of the second cell, given as ascending particle IDs.
   */
  ExpectedLists expectedNeighbors2;
};

/**
 * Fixture shared by the parameterized AoSFunctor, SoAFunctorSingle and SoAFunctorVerlet tests, which all run over the
 * same set of test cases.
 */
class NeighborIdentificationFunctorListGenTest : public NeighborIdentificationFunctorTest,
                                                 public ::testing::WithParamInterface<ListGenerationParams> {};

/**
 * Fixture for the parameterized SoAFunctorPair tests.
 */
class NeighborIdentificationFunctorSoAPairTest : public NeighborIdentificationFunctorTest,
                                                 public ::testing::WithParamInterface<PairListGenerationParams> {};

/**
 * Fixture for the parameterized full-flow tests, parameterized over the cross product of the newton3 modes and the
 * ownership variants.
 */
class NeighborIdentificationFunctorFullFlowTest
    : public NeighborIdentificationFunctorTest,
      public ::testing::WithParamInterface<std::tuple<N3Mode, OwnershipVariant>> {};

/**
 * Fixture shared by the thread safety tests (except Verlet).
 */
class NeighborIdentificationFunctorThreadSafeTest : public NeighborIdentificationFunctorTest,
                                                    public ::testing::WithParamInterface<N3Mode> {};