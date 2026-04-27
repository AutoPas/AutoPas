/**
 * @file LinkedCellsConstructorTest.cpp
 * @author Alexander Glas
 * @date 27.04.26
 */

#include "LinkedCellsConstructorTest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>
#include <vector>

#include "autopas/utils/ArrayUtils.h"

namespace {

using CellContentSnapshot = std::vector<std::vector<std::tuple<size_t, autopas::OwnershipState>>>;

template <class ParticleCell>
void expectCellBlocksEquivalent(const autopas::internal::CellBlock3D<ParticleCell> &fromCellSizeFactor,
                                const autopas::internal::CellBlock3D<ParticleCell> &fromCellsPerDim,
                                const std::vector<std::array<double, 3>> &positionsToCheck) {
  EXPECT_EQ(fromCellSizeFactor.getCellsPerDimensionWithHalo(), fromCellsPerDim.getCellsPerDimensionWithHalo());
  EXPECT_EQ(fromCellSizeFactor.getCellsPerDimensionWithoutHalo(), fromCellsPerDim.getCellsPerDimensionWithoutHalo());
  EXPECT_EQ(fromCellSizeFactor.getCellsPerInteractionLength(), fromCellsPerDim.getCellsPerInteractionLength());
  // test getCellsPerDimensionWithoutHalo:
  for (int d = 0; d < 3; ++d) {
    EXPECT_EQ(fromCellSizeFactor.getCellsPerDimensionWithHalo()[d],
              fromCellSizeFactor.getCellsPerDimensionWithoutHalo()[d] +
                  2 * fromCellSizeFactor.getCellsPerInteractionLength());
  }
  EXPECT_EQ(fromCellSizeFactor.getNumCells(), fromCellsPerDim.getNumCells());
  EXPECT_EQ(fromCellSizeFactor.getFirstOwnedCellIndex(), fromCellsPerDim.getFirstOwnedCellIndex());
  EXPECT_EQ(fromCellSizeFactor.getLastOwnedCellIndex(), fromCellsPerDim.getLastOwnedCellIndex());

  EXPECT_THAT(fromCellSizeFactor.getCellLength(),
              testing::Pointwise(testing::DoubleNear(1e-14), fromCellsPerDim.getCellLength()));
  EXPECT_THAT(fromCellSizeFactor.getHaloBoxMin(),
              testing::Pointwise(testing::DoubleNear(1e-14), fromCellsPerDim.getHaloBoxMin()));
  EXPECT_THAT(fromCellSizeFactor.getHaloBoxMax(),
              testing::Pointwise(testing::DoubleNear(1e-14), fromCellsPerDim.getHaloBoxMax()));

  for (const auto &pos : positionsToCheck) {
    EXPECT_EQ(fromCellSizeFactor.get3DIndexOfPosition(pos), fromCellsPerDim.get3DIndexOfPosition(pos))
        << "Position: " << autopas::utils::ArrayUtils::to_string(pos);
    EXPECT_EQ(fromCellSizeFactor.get1DIndexOfPosition(pos), fromCellsPerDim.get1DIndexOfPosition(pos))
        << "Position: " << autopas::utils::ArrayUtils::to_string(pos);
  }

  std::set<size_t> representativeCellIndices{0, fromCellSizeFactor.getFirstOwnedCellIndex(),
                                             fromCellSizeFactor.getLastOwnedCellIndex(),
                                             fromCellSizeFactor.getNumCells() - 1};
  for (const auto cellIndex : representativeCellIndices) {
    const auto [minA, maxA] = fromCellSizeFactor.getCellBoundingBox(cellIndex);
    const auto [minB, maxB] = fromCellsPerDim.getCellBoundingBox(cellIndex);
    EXPECT_THAT(minA, testing::Pointwise(testing::DoubleNear(1e-14), minB));
    EXPECT_THAT(maxA, testing::Pointwise(testing::DoubleNear(1e-14), maxB));
  }
}

template <class LinkedCellsType>
CellContentSnapshot getCellContentSnapshot(LinkedCellsType &linkedCells) {
  CellContentSnapshot snapshot;
  snapshot.reserve(linkedCells.getCells().size());

  for (auto &cell : linkedCells.getCells()) {
    std::vector<std::tuple<size_t, autopas::OwnershipState>> particleData;
    particleData.reserve(cell.size());
    for (const auto &particle : cell) {
      particleData.emplace_back(particle.getID(), particle.getOwnershipState());
    }
    std::sort(particleData.begin(), particleData.end(), [](const auto &a, const auto &b) {
      const auto [idA, stateA] = a;
      const auto [idB, stateB] = b;
      if (idA != idB) {
        return idA < idB;
      }
      return toInt64(stateA) < toInt64(stateB);
    });
    snapshot.push_back(std::move(particleData));
  }

  return snapshot;
}

}  // namespace

std::array<unsigned long, 3> LinkedCellsConstructorTest::cellsPerDimFromCellSizeFactor(
    const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, double interactionLength,
    double cellSizeFactor) {
  std::array<unsigned long, 3> cellsPerDim{};
  for (size_t d = 0; d < 3; ++d) {
    const auto boxLength = boxMax[d] - boxMin[d];
    cellsPerDim[d] =
        std::max(static_cast<unsigned long>(std::floor(boxLength / (interactionLength * cellSizeFactor))), 1ul);
  }
  return cellsPerDim;
}

std::vector<std::array<double, 3>> LinkedCellsConstructorTest::representativePositions(
    const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
  const auto center =
      std::array<double, 3>{(boxMin[0] + boxMax[0]) / 2., (boxMin[1] + boxMax[1]) / 2., (boxMin[2] + boxMax[2]) / 2.};
  return {
      center,
      {boxMin[0], boxMin[1], boxMin[2]},
      {std::nextafter(boxMin[0], std::numeric_limits<double>::infinity()), center[1], center[2]},
      {std::nextafter(boxMin[0], -std::numeric_limits<double>::infinity()), center[1], center[2]},
      {std::nextafter(boxMax[0], -std::numeric_limits<double>::infinity()), center[1], center[2]},
      {boxMax[0], boxMax[1], boxMax[2]},
      {std::nextafter(boxMax[0], std::numeric_limits<double>::infinity()), center[1], center[2]},
      {center[0], boxMin[1], boxMax[2]},
      {boxMax[0], center[1], boxMin[2]},
  };
}

std::vector<ParticleFP64> LinkedCellsConstructorTest::representativeOwnedParticles(
    const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) {
  const std::array<double, 3> v{0., 0., 0.};
  const auto dx = (boxMax[0] - boxMin[0]) / 10.;
  const auto dy = (boxMax[1] - boxMin[1]) / 10.;
  const auto dz = (boxMax[2] - boxMin[2]) / 10.;
  return {
      {{boxMin[0] + dx, boxMin[1] + dy, boxMin[2] + dz}, v, 1},
      {{boxMin[0] + 3. * dx, boxMin[1] + 5. * dy, boxMin[2] + 7. * dz}, v, 2},
      {{std::nextafter(boxMax[0], -std::numeric_limits<double>::infinity()), boxMin[1] + 2. * dy, boxMin[2] + 8. * dz},
       v,
       3},
  };
}

std::vector<ParticleFP64> LinkedCellsConstructorTest::representativeHaloParticles(const std::array<double, 3> &boxMin,
                                                                                  const std::array<double, 3> &boxMax) {
  const std::array<double, 3> v{0., 0., 0.};
  const auto cx = (boxMin[0] + boxMax[0]) / 2.;
  const auto cy = (boxMin[1] + boxMax[1]) / 2.;
  const auto cz = (boxMin[2] + boxMax[2]) / 2.;
  return {
      {{std::nextafter(boxMin[0], -std::numeric_limits<double>::infinity()), cy, cz},
       v,
       4,
       autopas::OwnershipState::halo},
      {{std::nextafter(boxMax[0], std::numeric_limits<double>::infinity()), cy, cz},
       v,
       5,
       autopas::OwnershipState::halo},
      {{cx, std::nextafter(boxMin[1], -std::numeric_limits<double>::infinity()), cz},
       v,
       6,
       autopas::OwnershipState::halo},
      {{cx, cy, std::nextafter(boxMax[2], std::numeric_limits<double>::infinity())},
       v,
       7,
       autopas::OwnershipState::halo},
  };
}

TEST_F(LinkedCellsConstructorTest, CellBlock3DConstructorsProvideEquivalentObjects) {
  struct CellBlockCase {
    std::array<double, 3> boxMin;
    std::array<double, 3> boxMax;
    double interactionLength;
    double cellSizeFactor;
  };

  const std::vector<CellBlockCase> cases{
      {{0., 0., 0.}, {10., 10., 10.}, 1., 1.0},     {{0., 0., 0.}, {10., 10., 10.}, 1., 0.5},
      {{0., 0., 0.}, {10., 10., 10.}, 1., 2.0},     {{2. / 3., 0., 0.}, {1., .125, .125}, 0.03, 1.0},
      {{0., 0., 0.}, {58.5, 58.5, 58.5}, 3.0, 1.0},
  };

  for (const auto &testCase : cases) {
    std::vector<FMCell> cellsFromCellSizeFactor;
    std::vector<FMCell> cellsFromCellsPerDim;

    const auto cellsPerDim = cellsPerDimFromCellSizeFactor(testCase.boxMin, testCase.boxMax, testCase.interactionLength,
                                                           testCase.cellSizeFactor);

    autopas::internal::CellBlock3D<FMCell> fromCellSizeFactor(cellsFromCellSizeFactor, testCase.boxMin, testCase.boxMax,
                                                              testCase.interactionLength, testCase.cellSizeFactor);
    autopas::internal::CellBlock3D<FMCell> fromCellsPerDim(cellsFromCellsPerDim, testCase.boxMin, testCase.boxMax,
                                                           testCase.interactionLength, cellsPerDim);

    expectCellBlocksEquivalent(fromCellSizeFactor, fromCellsPerDim,
                               representativePositions(testCase.boxMin, testCase.boxMax));
  }
}

TEST_F(LinkedCellsConstructorTest, LinkedCellsConstructorsProvideEquivalentObjects) {
  struct LinkedCellsCase {
    std::array<double, 3> boxMin;
    std::array<double, 3> boxMax;
    double cutoff;
    double skin;
    double cellSizeFactor;
  };

  const std::vector<LinkedCellsCase> cases{
      {{0., 0., 0.}, {10., 10., 10.}, 1.0, 0.0, 1.0},
      {{0., 0., 0.}, {10., 10., 10.}, 1.0, 0.5, 0.5},
      {{0., 0., 0.}, {10., 10., 10.}, 1.5, 0.5, 2.0},
      {{2. / 3., 0., 0.}, {1., .125, .125}, 0.02, 0.01, 1.0},
  };

  for (const auto &testCase : cases) {
    const auto interactionLength = testCase.cutoff + testCase.skin;
    const auto cellsPerDim =
        cellsPerDimFromCellSizeFactor(testCase.boxMin, testCase.boxMax, interactionLength, testCase.cellSizeFactor);

    autopas::LinkedCells<ParticleFP64> fromCellSizeFactor(testCase.boxMin, testCase.boxMax, testCase.cutoff,
                                                          testCase.skin, testCase.cellSizeFactor);
    autopas::LinkedCells<ParticleFP64> fromCellsPerDim(testCase.boxMin, testCase.boxMax, testCase.cutoff, testCase.skin,
                                                       cellsPerDim);

    const auto positions = representativePositions(testCase.boxMin, testCase.boxMax);
    expectCellBlocksEquivalent(fromCellSizeFactor.getCellBlock(), fromCellsPerDim.getCellBlock(), positions);

    for (const auto &particle : representativeOwnedParticles(testCase.boxMin, testCase.boxMax)) {
      fromCellSizeFactor.addParticle(particle);
      fromCellsPerDim.addParticle(particle);
    }
    for (const auto &particle : representativeHaloParticles(testCase.boxMin, testCase.boxMax)) {
      fromCellSizeFactor.addHaloParticle(particle);
      fromCellsPerDim.addHaloParticle(particle);
    }

    EXPECT_EQ(fromCellSizeFactor.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHaloOrDummy),
              fromCellsPerDim.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHaloOrDummy));

    EXPECT_EQ(getCellContentSnapshot(fromCellSizeFactor), getCellContentSnapshot(fromCellsPerDim));
  }
}