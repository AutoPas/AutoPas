/**
 * @file TraversalTest.cpp
 * @author C. Menges
 * @date 16.04.2019
 */

#include "TraversalTest.h"

#include "autopas/containers/CompatibleLoadEstimators.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/containers/cellPairTraversals/BalancedTraversal.h"
#include "autopas/selectors/TraversalSelector.h"
#include "autopas/selectors/TraversalSelectorInfo.h"
#include "testingHelpers/NumThreadGuard.h"

using ::testing::_;  // anything is ok
using ::testing::Bool;
using ::testing::Combine;
using ::testing::ValuesIn;

void testTraversal(autopas::TraversalOption traversalOption, autopas::LoadEstimatorOption loadEstimatorOption,
                   bool useN3, const std::array<size_t, 3> &edgeLength, int interactions, double cutoff = 1.0) {
  const std::array<double, 3> linkedCellsBoxMax = {(double)(edgeLength[0]), (double)(edgeLength[1]),
                                                   (double)(edgeLength[2])};
  const std::array<double, 3> linkedCellsBoxMin = {0., 0., 0.};

  TraversalTest::CountFunctor functor(cutoff);
  autopas::LinkedCells<FPCell> linkedCells(linkedCellsBoxMin, linkedCellsBoxMax, cutoff, 0.0, 1.0 / cutoff,
                                           loadEstimatorOption);

  autopasTools::generators::GridGenerator::fillWithParticles(linkedCells, edgeLength);
  ASSERT_EQ(linkedCells.getNumParticles(), edgeLength[0] * edgeLength[1] * edgeLength[2]);

  std::array<unsigned long, 3> overlap = {};
  for (unsigned int d = 0; d < 3; d++) {
    overlap[d] = std::ceil(cutoff / 1.0);
  }
  const auto cellsPerDim =
      autopas::utils::ArrayMath::add(edgeLength, autopas::utils::ArrayMath::mulScalar(overlap, 2ul));
  NumThreadGuard numThreadGuard(4);
  // clustersize is 32 if traversal has something like cluster in it, otherwise 0.
  unsigned int clusterSize = traversalOption.to_string().find("luster") != std::string::npos ? 32 : 0;
  // this test assumes a cell size of 1. in each direction
  autopas::TraversalSelectorInfo tsi(cellsPerDim, cutoff, {1., 1., 1.}, clusterSize);
  std::unique_ptr<autopas::TraversalInterface> traversal;
  if (useN3 and traversalOption != autopas::TraversalOption::lc_c01) {
    traversal = autopas::TraversalSelector<FPCell>::template generateTraversal<TraversalTest::CountFunctor,
                                                                               autopas::DataLayoutOption::aos, true>(
        traversalOption, functor, tsi);
  } else {
    traversal = autopas::TraversalSelector<FPCell>::template generateTraversal<TraversalTest::CountFunctor,
                                                                               autopas::DataLayoutOption::aos, false>(
        traversalOption, functor, tsi);
  }

  unsigned long cellId = 0;

  const auto boxMax = autopas::utils::ArrayMath::sub(edgeLength, overlap);

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

  auto *traversalInterface = traversal.get();
  linkedCells.iteratePairwise(traversalInterface);
}

TEST_P(TraversalTest, testTraversal_2x2x2) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {2ul, 2ul, 2ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_2x3x4) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {2ul, 3ul, 4ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_3x3x3) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {3ul, 3ul, 3ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_8x8x8) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {8ul, 8ul, 8ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_8x8x8_overlap2) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {8ul, 8ul, 8ul};
  const auto cutoff = 2.0;

  if (traversalOption == autopas::TraversalOption::lc_c04) {
    GTEST_SKIP_("C04 doesn't support cellSizeFactors < 1.0");
  }
  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 32, cutoff);
}

TEST_P(TraversalTest, testTraversal_6x7x8) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {6ul, 7ul, 8ul};
  const auto cutoff = 1.0;

  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 6, cutoff);
}

TEST_P(TraversalTest, testTraversal_6x7x8_overlap2) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {6ul, 7ul, 8ul};
  const auto cutoff = 2.0;

  if (traversalOption == autopas::TraversalOption::lc_c04) {
    GTEST_SKIP_("C04 doesn't support cellSizeFactors < 1.0");
  }
  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 32, cutoff);
}

TEST_P(TraversalTest, testTraversal_7x8x9_overlap3) {
  autopas::TraversalOption traversalOption = std::get<0>(GetParam());
  autopas::LoadEstimatorOption loadEstimatorOption = std::get<2>(GetParam());
  auto newton3 = std::get<1>(GetParam());
  std::array<size_t, 3> domain = {7ul, 8ul, 9ul};
  const auto cutoff = 3.0;

  if (traversalOption == autopas::TraversalOption::lc_c04) {
    GTEST_SKIP_("C04 doesn't support cellSizeFactors < 1.0");
  }
  testTraversal(traversalOption, loadEstimatorOption, newton3, domain, 122, cutoff);
}

INSTANTIATE_TEST_SUITE_P(
    Generated, TraversalTest,
    ValuesIn([]() -> std::set<std::tuple<autopas::TraversalOption, bool, autopas::LoadEstimatorOption>> {
      auto allTraversals = autopas::compatibleTraversals::allLCCompatibleTraversals();
      allTraversals.erase(autopas::TraversalOption::lc_c01_cuda);
      allTraversals.erase(autopas::TraversalOption::lc_c01_combined_SoA);
      allTraversals.erase(autopas::TraversalOption::lc_c04_combined_SoA);
      std::set<std::tuple<autopas::TraversalOption, bool, autopas::LoadEstimatorOption>> res;
      for (const auto &traversal : allTraversals) {
        for (const auto &loadEstimator : autopas::loadEstimators::getApplicableLoadEstimators(
                 autopas::ContainerOption::linkedCells, traversal, autopas::LoadEstimatorOption::getAllOptions())) {
          res.emplace(traversal, false, loadEstimator);
          res.emplace(traversal, true, loadEstimator);
        }
      }
      return res;
    }()),
    TraversalTest::PrintToStringParamName());
