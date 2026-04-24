/**
 * @file LJFunctorHWYBench.cpp
 *
 * Micro-benchmarks for every public entry point of mdLib::LJFunctorHWY:
 *   - AoSFunctor
 *   - SoAFunctorSingle
 *   - SoAFunctorPair (all four VectorizationPatterns)
 *   - SoAFunctorPairSorted
 *   - SoAFunctorVerlet
 *
 * Each case is parameterised by particle count and (where meaningful) by the
 * newton3 flag so both code paths are exercised. Particle data is generated
 * once per benchmark instance; force accumulators are reset between iterations
 * to keep the numeric range bounded.
 */

#include <benchmark/benchmark.h>

#include <array>
#include <random>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
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

/// Default template parameters used for every benchmark: single-site LJ,
/// shifted potential, no mixing, Newton3 supported, globals disabled so the
/// reduction does not pollute the kernel measurement, no FLOP counting.
using BenchFunctor = mdLib::LJFunctorHWY<Molecule, /*shifting=*/true, /*useMixing=*/false,
                                         autopas::FunctorN3Modes::Both, /*calculateGlobals=*/false,
                                         /*countFLOPs=*/false>;

/// Fill a single cell uniformly.
void fillCell(FMCell &cell, const std::array<double, 3> &low, const std::array<double, 3> &high, std::size_t n,
              unsigned seed) {
  const Molecule defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell, defaultParticle, low, high, n, seed);
}

/// Build a plain O(n^2) neighbor list within one cell for the Verlet path.
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

/// Reset force accumulators on an SoA buffer. Keeping forces bounded prevents
/// drift when the same kernel is run thousands of times on identical inputs.
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
    // reset forces each iteration
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
  // Report "pair evaluations per second" for a more intuitive metric.
  const double pairs = newton3 ? static_cast<double>(n) * (n - 1) / 2.0 : static_cast<double>(n) * (n - 1);
  state.counters["pairs/s"] = benchmark::Counter(pairs, benchmark::Counter::kIsIterationInvariantRate);
}
BENCHMARK(BM_AoSFunctor)
    ->ArgsProduct({{32, 128, 512}, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_AoS_Functor_Single");

// -----------------------------------------------------------------------------
// SoAFunctorSingle: one cell, self-interaction (internal newton3-like loop).
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
    ->ArgsProduct({{32, 128, 512, 2048}, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Functor_Single");

// -----------------------------------------------------------------------------
// SoAFunctorPair: two cells side by side, all four vectorization patterns.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPair(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  fillCell(cell1, kLow, {kHigh[0] / 2.0, kHigh[1], kHigh[2]}, n, 42);
  fillCell(cell2, {kHigh[0] / 2.0, kLow[1], kLow[2]}, kHigh, n, 1337);

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
    ->ArgsProduct({{32, 128, 512},
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Functor_Pair");

// -----------------------------------------------------------------------------
// SoAFunctorPairSorted: Gonnet-style pre-pruned SIMD. Only p1xVec is specialised
// internally; other patterns fall back to SoAFunctorPair.
// -----------------------------------------------------------------------------
static void BM_SoAFunctorPairSorted(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell1, cell2;
  fillCell(cell1, kLow, {kHigh[0] / 2.0, kHigh[1], kHigh[2]}, n, 42);
  fillCell(cell2, {kHigh[0] / 2.0, kLow[1], kLow[2]}, kHigh, n, 1337);

  auto functor = makeFunctor();
  functor.setVecPattern(VectorizationPattern::p1xVec);
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
    ->ArgsProduct({{32, 128, 512}, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Functor_SortedPair");

// -----------------------------------------------------------------------------
// SoAFunctorVerlet: full verlet sweep over an O(n^2) built neighbor list.
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
    ->ArgsProduct({{32, 128, 512}, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Functor_Verlet");

BENCHMARK_MAIN();
