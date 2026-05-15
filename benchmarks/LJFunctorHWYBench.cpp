/**
 * @file LJFunctorHWYBench.cpp
 *
 * Micro-benchmarks for mdLib::LJFunctorHWY, organised around two focused questions:
 *
 *   1. Ratio study — how does the interaction ratio (fraction of pairs within cutoff)
 *      affect SoAFunctorPair vs SoAFunctorPairSorted?
 *      VectorizationPattern fixed to p1xVec; ratio and N swept.
 *
 *   2. VecPattern study — how do the four VectorizationPatterns compare for
 *      SoAFunctorPair and SoAFunctorPairSorted?
 *      Ratio fixed at 50 %; vecPattern and N swept.
 *
 *   3. Context baselines — AoSFunctor, SoAFunctorSingle, SoAFunctorVerlet
 *      with varying N and newton3.
 *
 * Approximate run count: ~1050 (down from ~4730 in the previous layout).
 */

#include <benchmark/benchmark.h>

#include <array>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/TwoCellsInteractionRatioGenerator.h"
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
constexpr std::array<double, 3> kHigh{8.0, 8.0, 8.0};
constexpr double kBoundary = 8.0;   // x-coordinate of the shared cell boundary
constexpr double kCellSizeX = 8.0;  // x-extent of each pair cell (> 1.5 * kCutoff = 4.5)

// N values shared across all benchmarks: spans the interesting range without
// dense sampling where the scaling curve is already flat.
const std::vector<int64_t> kNValues = {10, 25, 50, 75, 100, 200, 250};

// vecPattern used for the ratio-study benchmarks (fixed so ratio is the only variable).
constexpr int kFixedVecPattern = static_cast<int>(VectorizationPattern::p1xVec);

// Interaction ratio (integer percent) used for the vecPattern-study benchmarks.
constexpr int kFixedRatioPct = 50;

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

template <class SoA>
void resetForces(SoA &soa) {
  auto *fx = soa.template begin<Molecule::AttributeNames::forceX>();
  auto *fy = soa.template begin<Molecule::AttributeNames::forceY>();
  auto *fz = soa.template begin<Molecule::AttributeNames::forceZ>();
  for (std::size_t i = 0; i < soa.size(); ++i) {
    fx[i] = 0.0;
    fy[i] = 0.0;
    fz[i] = 0.0;
  }
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
    for (auto &p : cell) p.setF({0.0, 0.0, 0.0});
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
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize=*/false);
  functor.initTraversal();

  for (auto _ : state) {
    resetForces(cell._particleSoABuffer);
    functor.SoAFunctorSingle(cell._particleSoABuffer, newton3);
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
// SoAFunctorPair — ratio study and vecPattern study share one benchmark body.
// Registration 1: vecPattern fixed, ratio swept  (core contribution comparison)
// Registration 2: ratio fixed at 50%, vecPattern swept  (pattern sensitivity)
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPair(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));
  const double ratio = state.range(3) / 100.0;

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::TwoCellsInteractionRatioGenerator::fillWithParticles(
      cell1, cell2, kLow, {kBoundary, kHigh[1], kHigh[2]}, {kBoundary, kLow[1], kLow[2]},
      {kBoundary + kCellSizeX, kHigh[1], kHigh[2]}, n, ratio, kCutoff, defaultParticle, /*seed=*/42);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
  functor.initTraversal();

  for (auto _ : state) {
    resetForces(cell1._particleSoABuffer);
    resetForces(cell2._particleSoABuffer);
    functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorPair)
    ->ArgsProduct({kNValues, {0, 1}, {kFixedVecPattern}, {0, 25, 50, 75, 100}})
    ->ArgNames({"N", "n3", "vecPat", "ratio%"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Ratio");
BENCHMARK(BM_SoAFunctorPair)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2), static_cast<int>(VectorizationPattern::pVecx1)},
                   {kFixedRatioPct}})
    ->ArgNames({"N", "n3", "vecPat", "ratio%"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_VecPatterns");

// -----------------------------------------------------------------------------
// SoAFunctorPairSorted — same two-registration pattern as SoAFunctorPair.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairSorted(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));
  const double ratio = state.range(3) / 100.0;

  FMCell cell1, cell2;
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::TwoCellsInteractionRatioGenerator::fillWithParticles(
      cell1, cell2, kLow, {kBoundary, kHigh[1], kHigh[2]}, {kBoundary, kLow[1], kLow[2]},
      {kBoundary + kCellSizeX, kHigh[1], kHigh[2]}, n, ratio, kCutoff, defaultParticle, /*seed=*/42);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
  functor.initTraversal();

  const std::array<double, 3> sortingDirection{1.0, 0.0, 0.0};

  for (auto _ : state) {
    resetForces(cell1._particleSoABuffer);
    resetForces(cell2._particleSoABuffer);
    functor.SoAFunctorPairSorted(cell1._particleSoABuffer, cell2._particleSoABuffer, sortingDirection, kCutoff,
                                 newton3);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
BENCHMARK(BM_SoAFunctorPairSorted)
    ->ArgsProduct({kNValues, {0, 1}, {kFixedVecPattern}, {0, 25, 50, 75, 100}})
    ->ArgNames({"N", "n3", "vecPat", "ratio%"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Ratio");
BENCHMARK(BM_SoAFunctorPairSorted)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2), static_cast<int>(VectorizationPattern::pVecx1)},
                   {kFixedRatioPct}})
    ->ArgNames({"N", "n3", "vecPat", "ratio%"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_VecPatterns");

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
  functor.SoALoader(cell, cell._particleSoABuffer, 0, false);
  functor.initTraversal();

  for (auto _ : state) {
    resetForces(cell._particleSoABuffer);
    for (std::size_t i = 0; i < n; ++i) {
      functor.SoAFunctorVerlet(cell._particleSoABuffer, i, lists[i], newton3);
    }
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
