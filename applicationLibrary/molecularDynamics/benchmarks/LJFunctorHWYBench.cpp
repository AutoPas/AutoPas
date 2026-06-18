/**
 * @file LJFunctorHWYBench.cpp
 *
 * Micro-benchmarks for mdLib::LJFunctorHWY covering three groups:
 *
 *   1. Face pair — hitrate-controlled (TwoCellsInteractionHitrateGenerator):
 *        a. Hitrate study: vecPattern fixed to p1xVec; hitrate and N swept.
 *        b. VecPattern study: all four vecPatterns × hitrate × N (reduced N).
 *      Cells are face-adjacent (shared yz-plane at x = kBoundary).
 *
 *   2. Geometry pair — uniform generator.
 *      Each geometry has an unsorted (SoAFunctorPair) and a sorted
 *      (SoAFunctorPairSorted) variant. all four vecPatterns are swept..
 *
 *   3. Baselines: AoSFunctor, SoAFunctorSingle, SoAFunctorVerlet — single cell,
 *      varying N and newton3. Reference points for absolute throughput.
 *
 * Cell size: kCellSize = 5 (≈ 1.67 × cutoff), N ∈ {10…150} per cell.
 */

#include <benchmark/benchmark.h>

#include <array>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/TwoCellsInteractionHitrateGenerator.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

namespace {

using Molecule = mdLib::MoleculeLJ;
using FMCell = autopas::FullParticleCell<Molecule>;
using VectorizationPattern = autopas::VectorizationPatternOption::Value;

constexpr double kCutoff = 3.0;
constexpr double kEpsilon = 1.0;
constexpr double kSigma = 1.0;
constexpr std::array<double, 3> kLow{0.0, 0.0, 0.0};
constexpr std::array<double, 3> kHigh{5.0, 5.0, 5.0};
constexpr double kBoundary = 5.0;  // x-coordinate of the shared cell boundary
constexpr double kCellSize = 5.0;  // x-extent of each pair cell (> 1.5 * kCutoff = 4.5)

const std::vector<int64_t> kNValues = {10, 25, 50, 75, 100, 125, 150};
const std::vector<int64_t> kNValuesReduced = {50, 75, 100};

const std::vector<int64_t> hitrates = {0, 5, 10, 15, 20, 30, 50};
const std::vector<int64_t> reducedHitrates = {5, 10, 15, 30};

// vecPattern used for the hitrate-study benchmarks (fixed so hitrate is the only variable).
constexpr int kFixedVecPattern = static_cast<int>(VectorizationPattern::p1xVec);

using BenchFunctor = mdLib::LJFunctorHWY<Molecule, /*shifting=*/true, /*useMixing=*/false,
                                         autopas::FunctorN3Modes::Both, /*calculateGlobals=*/false,
                                         /*countFLOPs=*/false>;

void fillCell(FMCell &cell, const std::array<double, 3> &low, const std::array<double, 3> &high, std::size_t n,
              unsigned seed) {
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell, defaultParticle, low, high, n, seed);
}

std::vector<std::vector<std::size_t, autopas::AlignedAllocator<std::size_t>>> buildNeighborLists(const FMCell &cell,
                                                                                                 double interactionLen,
                                                                                                 bool newton3) {
  using namespace autopas::utils::ArrayMath::literals;
  std::vector<std::vector<std::size_t, autopas::AlignedAllocator<std::size_t>>> lists(cell.size());
  const double interactionLenSq = interactionLen * interactionLen;
  for (std::size_t i = 0; i < cell.size(); ++i) {
    for (std::size_t j = newton3 ? i + 1 : 0; j < cell.size(); ++j) {
      if (i == j) continue;
      const auto dr = cell[i].getR() - cell[j].getR();
      if (autopas::utils::ArrayMath::dot(dr, dr) <= interactionLenSq) {
        lists[i].push_back(j);
      }
    }
  }
  return lists;
}

BenchFunctor makeFunctor() {
  BenchFunctor f(kCutoff);
  f.setParticleProperties(kEpsilon * 24.0, kSigma * kSigma);
  return f;
}

}  // namespace

// -----------------------------------------------------------------------------
// AoSFunctor: all O(n^2) pairs within one cell.
// -----------------------------------------------------------------------------
static void BM_AoSFunctor(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell;
  fillCell(cell, kLow, kHigh, n, 42);

  auto functor = makeFunctor();
  functor.initTraversal();

  for (auto _ : state) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = newton3 ? i + 1 : 0; j < n; ++j) {
        if (i == j) continue;
        functor.AoSFunctor(cell[i], cell[j], newton3);
      }
    }
    benchmark::DoNotOptimize(cell);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_AoSFunctor)->ArgsProduct({kNValues, {0, 1}})->ArgNames({"N", "n3"})->Repetitions(5)->Name("BM_AoS");

// -----------------------------------------------------------------------------
// SoAFunctorSingle: one cell, self-interaction.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorSingle(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell;
  fillCell(cell, kLow, kHigh, n, 42);

  auto functor = makeFunctor();
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize=*/false);
    functor.SoAFunctorSingle(cell._particleSoABuffer, newton3);
    functor.SoAExtractor(cell, cell._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorSingle)
    ->ArgsProduct({kNValues, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Single");

// -----------------------------------------------------------------------------
// Face pair, hitrate-controlled — SoAFunctorPair.
// Registration 1: vecPattern fixed (p1xVec), hitrate swept.
// Registration 2: all vecPatterns × hitrate (reduced N).
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairFace(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));
  const double hitrate = state.range(3) / 100.0;

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      cell1, cell2, kLow, {kBoundary, kHigh[1], kHigh[2]}, {kBoundary, kLow[1], kLow[2]},
      {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, hitrate, kCutoff, defaultParticle, /*seed=*/42);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorPairFace)
    ->ArgsProduct({kNValues, {0, 1}, {kFixedVecPattern}, hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Face_Hitrate");
BENCHMARK(BM_SoAFunctorPairFace)
    ->ArgsProduct({kNValuesReduced,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2), static_cast<int>(VectorizationPattern::pVecx1)},
                   hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Face_VecPatterns");

// -----------------------------------------------------------------------------
// Face pair, hitrate-controlled — SoAFunctorPairSorted. Same two registrations.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairSortedFace(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));
  const double hitrate = state.range(3) / 100.0;

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      cell1, cell2, kLow, {kBoundary, kHigh[1], kHigh[2]}, {kBoundary, kLow[1], kLow[2]},
      {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, hitrate, kCutoff, defaultParticle, /*seed=*/42);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  const std::array<double, 3> sortingDirection{1.0, 0.0, 0.0};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPairSorted(cell1._particleSoABuffer, cell2._particleSoABuffer, sortingDirection, kCutoff,
                                 newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorPairSortedFace)
    ->ArgsProduct({kNValues, {0, 1}, {kFixedVecPattern}, hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Face_Hitrate");
BENCHMARK(BM_SoAFunctorPairSortedFace)
    ->ArgsProduct({kNValuesReduced,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2), static_cast<int>(VectorizationPattern::pVecx1)},
                   hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Face_VecPatterns");

// -----------------------------------------------------------------------------
// Face pair, uniform generator — geometry-fixed hitrate (~30 %).
// Sweeps N and vecPattern only.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairFaceUniform(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, {kBoundary, kLow[1], kLow[2]},
                                                                {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorPairFaceUniform)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Face_Uniform");

static void BM_SoAFunctorSortedPairFaceUniform(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, {kBoundary, kLow[1], kLow[2]},
                                                                {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  const std::array<double, 3> sortingDirection{1.0, 0.0, 0.0};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPairSorted(cell1._particleSoABuffer, cell2._particleSoABuffer, sortingDirection, kCutoff,
                                 newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorSortedPairFaceUniform)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Face_Uniform");

// -----------------------------------------------------------------------------
// Edge pair, uniform generator — geometry-fixed hitrate (~10 %).
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairEdge(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2]},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2]}, n, 43);
  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

BENCHMARK(BM_SoAFunctorPairEdge)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Edge");

// -----------------------------------------------------------------------------
// Corner pair, uniform generator — geometry-fixed hitrate (~5 %).
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairCorner(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2] + kCellSize},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2] + kCellSize}, n, 43);
  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

BENCHMARK(BM_SoAFunctorPairCorner)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Corner");

static void BM_SoAFunctorSortedPairEdge(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2]},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  double normalized = 1.0 / sqrt(2.0);

  const std::array<double, 3> sortingDirection{normalized, normalized, 0.0};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPairSorted(cell1._particleSoABuffer, cell2._particleSoABuffer, sortingDirection, kCutoff,
                                 newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

BENCHMARK(BM_SoAFunctorSortedPairEdge)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Edge");

static void BM_SoAFunctorSortedPairCorner(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2] + kCellSize},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2] + kCellSize}, n, 43);
  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  const double invSqrt3 = 1.0 / sqrt(3.0);
  const std::array<double, 3> sortingDirection{invSqrt3, invSqrt3, invSqrt3};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    functor.SoAFunctorPairSorted(cell1._particleSoABuffer, cell2._particleSoABuffer, sortingDirection, kCutoff,
                                 newton3);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

BENCHMARK(BM_SoAFunctorSortedPairCorner)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Corner");

// -----------------------------------------------------------------------------
// SoAFunctorVerlet: Verlet sweep over an O(n^2)-built neighbor list.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorVerlet(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell;
  fillCell(cell, kLow, kHigh, n, 42);

  const double skin = 0.5;
  const auto lists = buildNeighborLists(cell, kCutoff + skin, newton3);

  auto functor = makeFunctor();
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell, cell._particleSoABuffer, 0, false);
    for (std::size_t i = 0; i < n; ++i) {
      functor.SoAFunctorVerlet(cell._particleSoABuffer, i, lists[i], newton3);
    }
    functor.SoAExtractor(cell, cell._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorVerlet)
    ->ArgsProduct({kNValues, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Verlet");

BENCHMARK_MAIN();
