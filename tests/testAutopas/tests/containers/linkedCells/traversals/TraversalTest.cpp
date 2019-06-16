/**
 * @file TraversalTest.cpp
 * @author C. Menges
 * @date 16.04.2019
 */

#include "TraversalTest.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/utils/Logger.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;  // anything is ok
using ::testing::Bool;
using ::testing::Combine;
using ::testing::ValuesIn;

void testTraversal(autopas::TraversalOption traversalOption, bool useN3, const std::array<size_t, 3> &edgeLength,
                   int interactions, double cutoff = 1.0) {
  TraversalTest::CountFunctor functor;
  functor.setCutoff(cutoff);
  std::vector<FPCell> cells(edgeLength[0] * edgeLength[1] * edgeLength[2]);

  GridGenerator::fillWithParticles<autopas::Particle>(cells, edgeLength);

  NumThreadGuard(4);

  autopas::TraversalSelectorInfo<FPCell> tsi(edgeLength, cutoff);
  std::unique_ptr<autopas::CellPairTraversal<FPCell>> Traversal;
  if (useN3 and traversalOption != autopas::TraversalOption::c01) {
    Traversal = autopas::TraversalSelector<FPCell>::template generateTraversal<TraversalTest::CountFunctor,
                                                                               autopas::DataLayoutOption::aos, true>(
        traversalOption, functor, tsi);
  } else {
    Traversal = autopas::TraversalSelector<FPCell>::template generateTraversal<TraversalTest::CountFunctor,
                                                                               autopas::DataLayoutOption::aos, false>(
        traversalOption, functor, tsi);
  }

  unsigned long cellId = 0;
  std::array<unsigned long, 3> overlap = {};
  for (unsigned int d = 0; d < 3; d++) {
    overlap[d] = std::ceil(cutoff / 1.0);
  }

  const auto boxMax = autopas::ArrayMath::sub(edgeLength, overlap);

  for (unsigned int z = 0; z < edgeLength[2]; ++z) {
    for (unsigned int y = 0; y < edgeLength[1]; ++y) {
      for (unsigned int x = 0; x < edgeLength[0]; ++x) {
        if (x >= overlap[0] && x < boxMax[0] && y >= overlap[1] && y < boxMax[1] && z >= overlap[2] && z < boxMax[2]) {
          EXPECT_CALL(functor, countFunc(cellId)).Times(interactions);
        } else {
          // Particles in the halo need to have less interactions as particles inside the box
          EXPECT_CALL(functor, countFunc(cellId)).Times(::testing::AtMost(interactions));
        }
        cellId++;
      }
    }
  }

  auto *traversalInterface = dynamic_cast<LinkedCellTraversalInterface<FPCell> *>(Traversal.get());
  traversalInterface->traverseCellPairs(cells);
}

TEST_P(TraversalTest, testTraversal_2x2x2) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {2ul, 2ul, 2ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_2x3x4) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {2ul, 3ul, 4ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_3x3x3) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {3ul, 3ul, 3ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_8x8x8) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {8ul, 8ul, 8ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_8x8x8_overlap2) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {8ul, 8ul, 8ul};
  const auto cutoff = 2.0;

  testTraversal(traversalOption, newton3, domain, 32, cutoff);
}

TEST_P(TraversalTest, testTraversal_6x7x8) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {6ul, 7ul, 8ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_6x7x8_overlap2) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {6ul, 7ul, 8ul};
  const auto cutoff = 2.0;

  testTraversal(traversalOption, newton3, domain, 32, cutoff);
}

TEST_P(TraversalTest, testTraversal_7x8x9_overlap3) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {7ul, 8ul, 9ul};
  const auto cutoff = 3.0;

  testTraversal(traversalOption, newton3, domain, 122, cutoff);
}

INSTANTIATE_TEST_SUITE_P(Generated, TraversalTest,
                         Combine(ValuesIn([]() -> std::set<autopas::TraversalOption> {
                                   auto allTraversals = autopas::compatibleTraversals::allLCCompatibleTraversals();
                                   allTraversals.erase(autopas::TraversalOption::c01Cuda);
                                   allTraversals.erase(autopas::TraversalOption::c01CombinedSoA);
                                   return allTraversals;
                                 }()),
                                 Bool()),
                         TraversalTest::PrintToStringParamName());
