/**
 * @file HGC08CellHandlerTest.cpp
 * @author Alexander Glas
 * @date 27.04.26
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "HGC08CellHandlerTest.h"
#include "autopas/containers/hierarchicalGrid/traversals/HGC08CellHandler.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::ElementsAre;
using ::testing::Pair;
using ::testing::UnorderedElementsAre;

namespace {

struct RecordingPairwiseFunctor : public EmptyPairwiseFunctor<Molecule> {
  std::vector<std::pair<size_t, size_t>> calls{};
  std::vector<bool> newton3Flags{};

  void AoSFunctor(Molecule &particle1, Molecule &particle2, bool newton3) override {
    calls.emplace_back(particle1.getID(), particle2.getID());
    newton3Flags.push_back(newton3);
  }
};

}  // namespace

void HGC08CellHandlerTest::runDecomposeScenario(const std::array<double, 3> &boxMin) {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> boxMax = boxMin + std::array<double, 3>{2., 2., 2.};

  std::vector<FMCell> lowerCells;
  std::vector<FMCell> upperCells;

  constexpr std::array<unsigned long, 3> lowerCellsPerDim{4, 4, 4};
  constexpr std::array<unsigned long, 3> upperCellsPerDim{2, 2, 2};
  const double cellLengthLower = 0.5;
  const double cellLengthUpper = 1.0;
  const double interactionLength = 0.5;

  autopas::internal::CellBlock3D<FMCell> lowerBlock(lowerCells, boxMin, boxMax, interactionLength, lowerCellsPerDim);
  autopas::internal::CellBlock3D<FMCell> upperBlock(upperCells, boxMin, boxMax, interactionLength, upperCellsPerDim);

  const std::array<size_t, 3> upperCell1Index{1, 1, 1};
  const std::array<size_t, 3> upperCell2Index{2, 1, 1};
  const auto upperCell1Center = (upperBlock.getCellBoundingBox(upperCell1Index).first +
                                 upperBlock.getCellBoundingBox(upperCell1Index).second) * 0.5;
  const auto upperCell2Center = (upperBlock.getCellBoundingBox(upperCell2Index).first +
                                 upperBlock.getCellBoundingBox(upperCell2Index).second) * 0.5;
  const auto upperCell1Index1D = upperBlock.get1DIndexOfPosition(upperCell1Center);
  const auto upperCell2Index1D = upperBlock.get1DIndexOfPosition(upperCell2Center);

  upperBlock.getCell(upperCell1Index).addParticle(Molecule(upperCell1Center, {0., 0., 0.}, 1));

  size_t particleId = 100;
  for (const size_t x : {3ul, 4ul}) {
    for (const size_t y : {1ul, 2ul}) {
      for (const size_t z : {1ul, 2ul}) {
        const auto index3D = std::array<size_t, 3>{x, y, z};
        const auto cellCenter = (lowerBlock.getCellBoundingBox(index3D).first + lowerBlock.getCellBoundingBox(index3D).second) * 0.5;
        lowerBlock.getCell(index3D).addParticle(Molecule(cellCenter, {0., 0., 0.}, particleId++));
      }
    }
  }

  RecordingPairwiseFunctor functor{};
  std::vector<autopas::internal::CellBlock3D<FMCell> *> cellBlocks{&lowerBlock, &upperBlock};
  autopas::HGC08CellHandler<FMCell, RecordingPairwiseFunctor> handler{
      &functor,
      upperCellsPerDim,
      interactionLength,
      std::array<double, 3>{cellLengthUpper, cellLengthUpper, cellLengthUpper},
      std::array<unsigned long, 3>{upperBlock.getCellsPerInteractionLength(), upperBlock.getCellsPerInteractionLength(),
                                   upperBlock.getCellsPerInteractionLength()},
      autopas::DataLayoutOption::aos,
      true,
      cellBlocks,
      std::vector<double>{0.2},
      1};

  handler.decompose2AndProcessCells(upperBlock.getCell(upperCell1Index), upperCell1Index1D, upperCell2Index1D);

  EXPECT_THAT(functor.calls,
              UnorderedElementsAre(Pair(1ul, 100ul), Pair(1ul, 101ul), Pair(1ul, 102ul), Pair(1ul, 103ul)));
  EXPECT_THAT(functor.newton3Flags, ElementsAre(true, true, true, true));
}

void HGC08CellHandlerTest::runPrecomputedOffsetsScenario() {
  using namespace autopas::utils::ArrayMath::literals;

  const std::array<double, 3> boxMin{0., 0., 0.};
  const std::array<double, 3> boxMax = boxMin + std::array<double, 3>{2., 2., 2.};

  std::vector<FMCell> lowerCells;
  std::vector<FMCell> upperCells;

  constexpr std::array<unsigned long, 3> lowerCellsPerDim{4, 4, 4};
  constexpr std::array<unsigned long, 3> upperCellsPerDim{2, 2, 2};
  const double cellLengthLower = 0.5;
  const double cellLengthUpper = 1.0;
  const double interactionLength = 0.5;

  autopas::internal::CellBlock3D<FMCell> lowerBlock(lowerCells, boxMin, boxMax, interactionLength, lowerCellsPerDim);
  autopas::internal::CellBlock3D<FMCell> upperBlock(upperCells, boxMin, boxMax, interactionLength, upperCellsPerDim);

  RecordingPairwiseFunctor functor{};
  std::vector<autopas::internal::CellBlock3D<FMCell> *> cellBlocks{&lowerBlock, &upperBlock};
  autopas::HGC08CellHandler<FMCell, RecordingPairwiseFunctor> handler{
      &functor,
      upperCellsPerDim,
      interactionLength,
      std::array<double, 3>{cellLengthUpper, cellLengthUpper, cellLengthUpper},
      std::array<unsigned long, 3>{upperBlock.getCellsPerInteractionLength(), upperBlock.getCellsPerInteractionLength(),
                                   upperBlock.getCellsPerInteractionLength()},
      autopas::DataLayoutOption::aos,
      true,
      cellBlocks,
      std::vector<double>{0.2},
      1};

  // Build reference set by re-running the decomposition logic for lowerLevel = 0
  std::vector<std::pair<unsigned long, unsigned long>> expectedPairs;
  const auto offsetsVec = autopas::LCC08CellHandlerUtility::computePairwiseCellOffsetsC08<
      autopas::LCC08CellHandlerUtility::C08OffsetMode::c08CellPairsSorting>(upperCellsPerDim,
                                                                           std::array<double, 3>{cellLengthUpper,
                                                                                                cellLengthUpper,
                                                                                                cellLengthUpper},
                                                                           interactionLength);
  for (auto const &tp : offsetsVec) {
    const auto offset1 = std::get<0>(tp);
    const auto offset2 = std::get<1>(tp);

    auto [lowCornerCell2, highCornerCell2] = upperBlock.getCellBoundingBox(offset2);
    // apply small shift like production code
    const std::array<double, 3> shift = lowerBlock.getCellLength() * 0.01;
    lowCornerCell2 = {lowCornerCell2[0] + shift[0], lowCornerCell2[1] + shift[1], lowCornerCell2[2] + shift[2]};
    highCornerCell2 = {highCornerCell2[0] - shift[0], highCornerCell2[1] - shift[1], highCornerCell2[2] - shift[2]};

    auto startIndex3D = lowerBlock.get3DIndexOfPosition(lowCornerCell2);
    auto stopIndex3D = lowerBlock.get3DIndexOfPosition(highCornerCell2);
    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          auto [lowCornerLowerCell, highCornerLowerCell] = lowerBlock.getCellBoundingBox({x, y, z});
          using autopas::utils::ArrayMath::boxDistanceSquared;
          const auto upperBBox1 = upperBlock.getCellBoundingBox(offset1);
          if (boxDistanceSquared(lowCornerLowerCell, highCornerLowerCell, upperBBox1.first, upperBBox1.second) <=
              interactionLength * interactionLength) {
            const auto lowerIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, lowerBlock.getCellsPerDimensionWithHalo());
            expectedPairs.emplace_back(offset1, lowerIndex);
          }
        }
      }
    }
  }

  // Extract actual pairs from handler
  std::vector<std::pair<unsigned long, unsigned long>> actualPairs;
  for (const auto &t : handler._cellPairOffsetsPerLevel[0]) {
    actualPairs.emplace_back(std::get<0>(t), std::get<1>(t));
  }

  std::sort(expectedPairs.begin(), expectedPairs.end());
  std::sort(actualPairs.begin(), actualPairs.end());
  EXPECT_EQ(expectedPairs, actualPairs);
}

TEST_F(HGC08CellHandlerTest, Decompose2AndProcessCellsSelectsExpectedLowerCells_ZeroOrigin) {
  runDecomposeScenario({0., 0., 0.});
}

TEST_F(HGC08CellHandlerTest, Decompose2AndProcessCellsSelectsExpectedLowerCells_ShiftedOrigin) {
  runDecomposeScenario({1., 2., 3.});
}

TEST_F(HGC08CellHandlerTest, PrecomputedInterLevelOffsetsMatchReference) {
  runPrecomputedOffsetsScenario();
}
