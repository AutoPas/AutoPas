/**
 * @file HGFittedTest.cpp
 * @author Alexander Glas
 * @date 27.04.26
 */

#include "HGFittedTest.h"

#include <gtest/gtest.h>

#include <array>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "autopas/containers/hierarchicalGrid/HierarchicalGrid.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "testingHelpers/commonTypedefs.h"

namespace {

struct ParsedLevel {
  std::array<double, 3> cellLength{};
  std::size_t cellsPerInteractionLength{};
  std::array<std::size_t, 3> cellsPerDimensionWithHalo{};
};

std::string extractAfterPrefix(const std::string &line, const std::string &prefix) {
  if (line.rfind(prefix, 0) != 0) {
    throw std::runtime_error("Unexpected line prefix: " + line);
  }
  return line.substr(prefix.size());
}

std::string extractBetweenMarkers(const std::string &line, const std::string &prefix, const std::string &suffix) {
  if (line.rfind(prefix, 0) != 0) {
    throw std::runtime_error("Unexpected line prefix: " + line);
  }

  const auto suffixPosition = line.find(suffix, prefix.size());
  if (suffixPosition == std::string::npos) {
    throw std::runtime_error("Unexpected line format: " + line);
  }

  return line.substr(prefix.size(), suffixPosition - prefix.size());
}

template <class T, std::size_t N>
std::array<T, N> parseArray(const std::string &text) {
  std::array<T, N> values{};
  std::istringstream stream(text);

  char openingBracket = '\0';
  if (!(stream >> openingBracket) || openingBracket != '[') {
    throw std::runtime_error("Expected array opening bracket: " + text);
  }

  for (std::size_t index = 0; index < N; ++index) {
    if (!(stream >> values[index])) {
      throw std::runtime_error("Could not parse array entry: " + text);
    }
    if (index + 1 < N) {
      char separator = '\0';
      if (!(stream >> separator) || separator != ',') {
        throw std::runtime_error("Expected array separator: " + text);
      }
    }
  }

  char closingBracket = '\0';
  if (!(stream >> closingBracket) || closingBracket != ']') {
    throw std::runtime_error("Expected array closing bracket: " + text);
  }

  return values;
}

ParsedLevel parseLevel(const std::string &cellLine, const std::string &interactionLine, const std::string &haloLine) {
  const std::string cellLengthPrefix{"CellLength: "};
  const std::string ciplPrefix{"CellsPerInteractionLength: "};
  const std::string numCellsMarker{" NumCells: "};
  const std::string haloPrefix{"CellsPerDimensionWithHalo: "};

  if (cellLine.rfind(cellLengthPrefix, 0) != 0) {
    throw std::runtime_error("Unexpected cell length line: " + cellLine);
  }

  ParsedLevel level{};
  level.cellLength = parseArray<double, 3>(extractAfterPrefix(cellLine, cellLengthPrefix));

  if (interactionLine.rfind(ciplPrefix, 0) != 0 || interactionLine.find(numCellsMarker) == std::string::npos) {
    throw std::runtime_error("Unexpected cells-per-interaction-length line: " + interactionLine);
  }

  level.cellsPerInteractionLength = std::stoul(extractBetweenMarkers(interactionLine, ciplPrefix, numCellsMarker));

  if (haloLine.rfind(haloPrefix, 0) != 0) {
    throw std::runtime_error("Unexpected halo line: " + haloLine);
  }
  level.cellsPerDimensionWithHalo = parseArray<std::size_t, 3>(extractAfterPrefix(haloLine, haloPrefix));

  return level;
}

std::array<double, 3> haloWidth(const ParsedLevel &level) {
  return {level.cellLength[0] * static_cast<double>(level.cellsPerInteractionLength),
          level.cellLength[1] * static_cast<double>(level.cellsPerInteractionLength),
          level.cellLength[2] * static_cast<double>(level.cellsPerInteractionLength)};
}

std::vector<ParsedLevel> parseLevels(const std::string &hierarchyDump) {
  std::istringstream stream(hierarchyDump);
  std::string line;
  std::vector<ParsedLevel> levels;

  while (std::getline(stream, line)) {
    if (line.rfind("HierarchicalGrid ", 0) != 0 || line.find(" numParticles: ") == std::string::npos) {
      continue;
    }

    std::string cellLine;
    std::string interactionLine;
    std::string haloLine;

    if (!std::getline(stream, cellLine) || !std::getline(stream, interactionLine) || !std::getline(stream, haloLine)) {
      throw std::runtime_error("Unexpected end of hierarchy dump while parsing level information.");
    }

    levels.push_back(parseLevel(cellLine, interactionLine, haloLine));
  }

  return levels;
}

std::size_t cellsPerDimensionWithoutHalo(const ParsedLevel &level, std::size_t dim) {
  return level.cellsPerDimensionWithHalo[dim] - 2 * level.cellsPerInteractionLength;
}

}  // namespace

void HGFittedTest::expectFittedHierarchy(const std::string &hierarchyDump, const std::array<double, 3> &boxMin,
                                         const std::array<double, 3> &boxMax) {
  const auto levels = parseLevels(hierarchyDump);
  ASSERT_GE(levels.size(), 2u) << "Expected at least two fitted hierarchy levels.";

  const std::array<double, 3> boxLength{boxMax[0] - boxMin[0], boxMax[1] - boxMin[1], boxMax[2] - boxMin[2]};
  constexpr double tolerance = 1e-5;

  for (std::size_t levelIndex = 0; levelIndex < levels.size(); ++levelIndex) {
    SCOPED_TRACE(::testing::Message() << "level " << levelIndex);

    const auto &level = levels[levelIndex];
    const auto width = haloWidth(level);

    for (std::size_t dim = 0; dim < 3; ++dim) {
      const auto cellsWithoutHalo = cellsPerDimensionWithoutHalo(level, dim);
      ASSERT_GT(cellsWithoutHalo, 0u);
      EXPECT_NEAR(level.cellLength[dim] * static_cast<double>(cellsWithoutHalo), boxLength[dim], tolerance);
      EXPECT_EQ(level.cellsPerDimensionWithHalo[dim], cellsWithoutHalo + 2 * level.cellsPerInteractionLength);
      EXPECT_GT(level.cellsPerInteractionLength, 0u);
      EXPECT_GT(width[dim], 0.0);
    }
  }

  for (std::size_t levelIndex = 0; levelIndex + 1 < levels.size(); ++levelIndex) {
    SCOPED_TRACE(::testing::Message() << "levels " << levelIndex << " and " << levelIndex + 1);

    const auto &finer = levels[levelIndex];
    const auto &coarser = levels[levelIndex + 1];
    const auto finerHaloWidth = haloWidth(finer);
    const auto coarserHaloWidth = haloWidth(coarser);

    for (std::size_t dim = 0; dim < 3; ++dim) {
      const auto finerCellsWithoutHalo = cellsPerDimensionWithoutHalo(finer, dim);
      const auto coarserCellsWithoutHalo = cellsPerDimensionWithoutHalo(coarser, dim);

      ASSERT_GT(coarserCellsWithoutHalo, 0u);
      EXPECT_EQ(finerCellsWithoutHalo % coarserCellsWithoutHalo, 0u)
          << "The refined grid does not fit an integer number of times into the coarser grid in dimension " << dim;
      EXPECT_EQ(finerCellsWithoutHalo * coarser.cellsPerInteractionLength,
                coarserCellsWithoutHalo * finer.cellsPerInteractionLength)
          << "The halo thickness ratio does not match the grid ratio in dimension " << dim;
      EXPECT_NEAR(finerHaloWidth[dim], coarserHaloWidth[dim], tolerance)
          << "The halo thickness differs between adjacent levels in dimension " << dim;
    }
  }
}

TEST_F(HGFittedTest, HierarchyLevelsFitIntoNextOneExactly) {
  struct TestCase {
    std::array<double, 3> boxMin;
    std::array<double, 3> boxMax;
    std::vector<double> cutoffs;
    double skin;
    double cellSizeFactor;
  };

  const std::vector<TestCase> testCases{
      {{0., 0., 0.}, {60., 60., 60.}, {0.5, 1.0, 2.0, 4.0}, 0.25, 1.0},
      {{-3., 2., 5.}, {33., 50., 65.}, {0.75, 1.5, 3.0}, 0.5, 2.0},
      {{1., 2., 3.}, {25., 32., 41.}, {1.0, 1.4}, 0.1, 0.5},
  };

  for (std::size_t caseIndex = 0; caseIndex < testCases.size(); ++caseIndex) {
    SCOPED_TRACE(::testing::Message() << "case " << caseIndex);

    const auto &testCase = testCases[caseIndex];
    autopas::HierarchicalGrid<ParticleFP64> container(
        testCase.boxMin, testCase.boxMax, testCase.cutoffs, testCase.skin, testCase.cellSizeFactor,
        /*sortingThreshold=*/8, autopas::LoadEstimatorOption::squaredParticlesPerCell,
        /*fittedGrids=*/true);

    const auto levels = parseLevels(container.toString());
    ASSERT_GE(levels.size(), 2u) << "Expected at least two fitted hierarchy levels.";

    for (std::size_t levelIndex = 0; levelIndex + 1 < levels.size(); ++levelIndex) {
      SCOPED_TRACE(::testing::Message() << "levels " << levelIndex << " and " << levelIndex + 1);

      const auto &finer = levels[levelIndex];
      const auto &coarser = levels[levelIndex + 1];

      for (std::size_t dim = 0; dim < 3; ++dim) {
        const auto finerCellsWithoutHalo = cellsPerDimensionWithoutHalo(finer, dim);
        const auto coarserCellsWithoutHalo = cellsPerDimensionWithoutHalo(coarser, dim);

        ASSERT_GT(coarserCellsWithoutHalo, 0u);
        EXPECT_EQ(finerCellsWithoutHalo % coarserCellsWithoutHalo, 0u)
            << "The refined grid does not fit an integer number of times into the coarser grid in dimension " << dim;
      }
    }
  }
}