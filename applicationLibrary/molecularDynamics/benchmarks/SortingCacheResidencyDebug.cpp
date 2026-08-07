/**
 * @file SortingCacheResidencyDebug.cpp
 * @date 2026-08-06
 *
 * Standalone debug tool for the sorting-threshold calibration investigation. Does NOT touch any
 * production code, it only re-implements the measurement loop of
 * src/autopas/utils/SortingThresholdBenchmark.h so that the loop's cache behaviour can be varied.
 *
 * Hypothesis under test: SortingThresholdBenchmark::executeRun() times _iterations back-to-back
 * calls on the SAME two cells, so after the first call both cells and the CellFunctor's sorted
 * buffers are cache resident. A real traversal touches each cell pair exactly once. The sorted
 * path's overhead is dominated by the extra SortedSoAView fill and the scatter-back, which is
 * exactly the cost the reuse hides, so the calibration should read a lower crossover than the
 * traversal actually needs.
 *
 * Two modes, identical in every other respect (same geometry, same scatter, same alternation, same
 * particle counts, same number of timed calls, same decision rule):
 *   hot  reproduces the production benchmark, kIterations calls on one reused cell pair.
 *   cold draws each timed call from a pool of distinct cell pairs of a given total footprint, so a
 *        call's data has been evicted before that pair comes round again.
 *
 * The pool is sized by BYTES, not by pair count, because a pair's footprint scales with the
 * particle count under test. A fixed pair count would leave the pool cache resident at small N and
 * the cold measurement would silently degrade into a second hot measurement there.
 *
 * Pool footprints are absolute and identical on every machine. How far up the footprint axis the
 * sorted path stays behind is a property of the machine's cache and is part of what is being
 * measured, so it is deliberately not normalised away by detecting cache size and scaling to it.
 *
 * Three programmes, see usage() below. Not part of the thesis-cited benchmark dataset. Debug branch
 * only.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/SortingThresholdInfoSingle.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

namespace {

using MoleculeType = mdLib::MoleculeLJ;
using BenchCell = autopas::FullParticleCell<MoleculeType>;
using BenchFunctor = mdLib::LJFunctorHWY<MoleculeType, /*shifting=*/true, /*useMixing=*/false,
                                         autopas::FunctorN3Modes::Both, /*calculateGlobals=*/false,
                                         /*countFLOPs=*/false>;
using BenchCF = autopas::internal::CellFunctor<BenchCell, BenchFunctor, /*bidirectional=*/true>;

/// Mirrors SortingThresholdBenchmark::_iterations.
constexpr size_t kIterations = 25;
/// Mirrors SortingThresholdBenchmark::_repetitions.
constexpr size_t kRepetitions = 100;
/// Mirrors SortingThresholdBenchmark::_sortedWinMarginFraction.
constexpr double kSortedWinMargin = 0.05;
/// Mirrors SortingThresholdBenchmark::_requiredSortedWinRatio.
constexpr double kRequiredWinRatio = 0.7;
/// Mirrors SortingThresholdBenchmark::_scatterFactor.
constexpr double kScatterFactor = 0.2;
/// Cutoff. The benchmark sizes its cells by exactly this, so the cell-edge to sorting-cutoff ratio is 1.
constexpr double kCutoff = 2.5;

/// Base seed for the split draws. Fixed so a run is reproducible, overridable to resample.
size_t gSeed = 0x5eed;
/// Cold-pool footprint in MiB for the full and itersweep programmes. Absolute, not cache-relative.
double gPoolMiB = 256.0;

/**
 * Bytes a single particle occupies across the AoS storage and the SoA buffer it is loaded into.
 * The SoA holds one vector per entry of SoAArraysType, each 8 B wide for this particle type. Vector
 * over-allocation and the CellFunctor's sorted scratch buffers are not counted, so this is a lower
 * bound and the realised footprint is somewhat larger than the requested one.
 */
constexpr size_t kSoAArrays = std::tuple_size_v<typename MoleculeType::SoAArraysType>;
constexpr size_t kBytesPerParticle = sizeof(MoleculeType) + kSoAArrays * sizeof(double);

struct Geometry {
  std::array<double, 3> cell1Low{0., 0., 0.};
  std::array<double, 3> cell1High{kCutoff, kCutoff, kCutoff};
  std::array<double, 3> cell2Low{};
  std::array<double, 3> cell2High{};
  std::array<double, 3> sortingDirection{};
};

/// Reproduces the cell placement of SortingThresholdBenchmark::executeRun().
Geometry makeGeometry(const std::string &direction) {
  const double invSqrt3 = 1. / std::sqrt(3.);
  const double invSqrt2 = 1. / std::sqrt(2.);
  Geometry g;
  if (direction == "corner") {
    g.cell2Low = {kCutoff, kCutoff, kCutoff};
    g.cell2High = {2. * kCutoff, 2. * kCutoff, 2. * kCutoff};
    g.sortingDirection = {invSqrt3, invSqrt3, invSqrt3};
  } else if (direction == "edge") {
    g.cell2Low = {kCutoff, kCutoff, 0.};
    g.cell2High = {2. * kCutoff, 2. * kCutoff, kCutoff};
    g.sortingDirection = {invSqrt2, invSqrt2, 0.};
  } else {
    g.cell2Low = {kCutoff, 0., 0.};
    g.cell2High = {2. * kCutoff, kCutoff, kCutoff};
    g.sortingDirection = {1., 0., 0.};
  }
  return g;
}

/// Splits a combined count into the two per-cell counts, reproducing the benchmark's scatter.
struct Split {
  size_t n1;
  size_t n2;
};

Split makeSplit(size_t numParticles, std::mt19937 &gen) {
  const auto numParticlesEqDistr = static_cast<size_t>(std::ceil(numParticles * (1. - kScatterFactor)));
  const size_t numParticlesScatter = numParticles - numParticlesEqDistr + numParticlesEqDistr % 2;
  const size_t numParticlesPerCell = numParticlesEqDistr / 2;
  std::uniform_int_distribution<size_t> distrib(0, numParticlesScatter);
  const size_t toAddCell1 = distrib(gen);
  return {numParticlesPerCell + toAddCell1, numParticlesPerCell + numParticlesScatter - toAddCell1};
}

/// One (cell1, cell2) sample, kept alive so the cold pool can hold many of them at once.
struct CellPair {
  BenchCell cell1;
  BenchCell cell2;
};

void fillPair(CellPair &pair, const Geometry &g, const MoleculeType &defaultParticle, const Split &split,
              unsigned int seed) {
  pair.cell1.clear();
  pair.cell2.clear();
  autopas::generators::UniformGenerator::fillWithParticles(pair.cell1, defaultParticle, g.cell1Low, g.cell1High,
                                                           split.n1, seed);
  autopas::generators::UniformGenerator::fillWithParticles(pair.cell2, defaultParticle, g.cell2Low, g.cell2High,
                                                           split.n2, seed + 1u);
}

/**
 * Number of pool pairs needed to reach a byte footprint at a given combined particle count.
 * Rounded up to a multiple of 2 * iterations so the pool's two halves stay aligned to the repeating
 * split pattern, and floored so each half can always supply one full timed region.
 */
size_t poolPairsForFootprint(size_t targetBytes, size_t numParticles, size_t iterations) {
  const size_t bytesPerPair = std::max<size_t>(1, numParticles * kBytesPerParticle);
  const size_t alignment = 2 * iterations;
  size_t pairs = (targetBytes + bytesPerPair - 1) / bytesPerPair;
  pairs = std::max(pairs, alignment);
  return ((pairs + alignment - 1) / alignment) * alignment;
}

/**
 * Result of one measurement point. Times are per timed call so points taken with different loop
 * splits stay comparable.
 */
struct RunResult {
  double nsPerCallSorted;
  double nsPerCallUnsorted;
  double winFraction;
  size_t poolPairs;
  double poolMiB;
};

/**
 * Measures the sorted and unsorted paths at one combined particle count.
 *
 * hot == true reproduces SortingThresholdBenchmark::executeRun(), including the reuse of one cell
 * pair across the timed region, the reuse of its allocation across repetitions, and the alternation
 * of which path runs first. hot == false performs the same number of timed calls against a pool of
 * distinct pairs sized to targetBytes, so no call reuses the previous call's data.
 *
 * In cold mode the pool's split pattern repeats with period iterations, and the sorted path starts
 * half a pool away from the unsorted one. Both paths therefore see an identical sequence of
 * particle counts while reading disjoint memory, which is what keeps the ratio meaningful.
 *
 * @param pool Scratch pool, grown as needed and reused across calls. Contents are overwritten.
 */
RunResult measure(BenchFunctor &functor, const MoleculeType &defaultParticle, const Geometry &g, size_t numParticles,
                  bool newton3, bool hot, std::vector<CellPair> &pool, size_t targetBytes, size_t seed,
                  size_t iterations = kIterations, size_t repetitions = kRepetitions) {
  BenchCF cellFunctor{functor, kCutoff, autopas::DataLayoutOption::soa, newton3};
  const autopas::SortingThresholdInfoSingle zeroThreshold(0);
  cellFunctor.setSoASortingThresholds(zeroThreshold);
  cellFunctor.setAoSSortingThresholds(zeroThreshold);

  autopas::utils::Timer sortedTimer, unsortedTimer;
  size_t sortedWins = 0;
  std::mt19937 gen(static_cast<std::mt19937::result_type>(seed));

  const size_t poolPairs = hot ? 0 : poolPairsForFootprint(targetBytes, numParticles, iterations);
  const size_t half = poolPairs / 2;
  size_t repOffset = 0;

  if (not hot) {
    if (pool.size() < poolPairs) {
      pool.resize(poolPairs);
    }
    std::vector<Split> splitTable(iterations);
    for (auto &s : splitTable) {
      s = makeSplit(numParticles, gen);
    }
    for (size_t p = 0; p < poolPairs; ++p) {
      fillPair(pool[p], g, defaultParticle, splitTable[p % iterations], static_cast<unsigned int>(2 * p));
    }
  }

  // Production declares its two cells once and clear()s them per repetition, so the allocation is
  // reused across repetitions. Matched here rather than constructing a fresh pair inside the loop.
  CellPair hotPair;

  for (size_t i = 0; i < repetitions; ++i) {
    if (hot) {
      fillPair(hotPair, g, defaultParticle, makeSplit(numParticles, gen), static_cast<unsigned int>(2 * i));
    }

    const auto runOne = [&](bool sorted, autopas::utils::Timer &timer) {
      const std::array<double, 3> dir = sorted ? g.sortingDirection : std::array<double, 3>{0., 0., 0.};
      const size_t base = hot ? 0 : (repOffset + (sorted ? half : 0)) % poolPairs;
      timer.start();
      for (size_t j = 0; j < iterations; ++j) {
        CellPair &pair = hot ? hotPair : pool[(base + j) % poolPairs];
        functor.SoALoader(pair.cell1, pair.cell1._particleSoABuffer, 0, false);
        functor.SoALoader(pair.cell2, pair.cell2._particleSoABuffer, 0, false);
        cellFunctor.processCellPair(pair.cell1, pair.cell2, dir);
      }
      return timer.stop();
    };

    long unsortedDelta = 0;
    long sortedDelta = 0;
    if (i % 2 == 0) {
      unsortedDelta = runOne(false, unsortedTimer);
      sortedDelta = runOne(true, sortedTimer);
    } else {
      sortedDelta = runOne(true, sortedTimer);
      unsortedDelta = runOne(false, unsortedTimer);
    }

    if (static_cast<double>(sortedDelta) < static_cast<double>(unsortedDelta) * (1. - kSortedWinMargin)) {
      ++sortedWins;
    }
    // Advance to pairs neither path has touched yet. Both halves advance together, so the workload
    // stays matched between the paths.
    if (not hot) {
      repOffset = (repOffset + iterations) % half;
    }
  }

  const auto calls = static_cast<double>(repetitions * iterations);
  const double poolMiB =
      hot ? 0.0 : static_cast<double>(poolPairs * numParticles * kBytesPerParticle) / (1024.0 * 1024.0);
  return {static_cast<double>(sortedTimer.getTotalTime()) / calls,
          static_cast<double>(unsortedTimer.getTotalTime()) / calls,
          static_cast<double>(sortedWins) / static_cast<double>(repetitions), poolPairs, poolMiB};
}

void usage() {
  std::printf(
      "usage: SortingCacheResidencyDebug [full|itersweep|poolsweep] [poolMiB] [seed]\n"
      "  full      sweep N in hot and cold mode, report the crossover under the production rule\n"
      "  itersweep hold total timed calls fixed and vary the loop split, against a cold reference\n"
      "  poolsweep vary the cold pool footprint to trace the reuse-distance curve\n"
      "  poolMiB   cold pool footprint for full and itersweep, ignored by poolsweep\n");
}

}  // namespace

int main(int argc, char **argv) {
  std::string mode = "full";
  if (argc > 1) {
    mode = argv[1];
  }
  if (mode == "-h" or mode == "--help") {
    usage();
    return 0;
  }
  if (argc > 2) {
    gPoolMiB = std::strtod(argv[2], nullptr);
  }
  if (argc > 3) {
    gSeed = std::strtoull(argv[3], nullptr, 10);
  }
  const auto targetBytes = static_cast<size_t>(gPoolMiB * 1024.0 * 1024.0);

  BenchFunctor functor(kCutoff);
  functor.setParticleProperties(24.0, 1.0);
  functor.initTraversal();
  const MoleculeType defaultParticle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

  const std::vector<std::string> directions{"corner", "edge", "face"};

  std::printf("# mode=%s cutoff=%.3f cellEdgeOverSortingCutoff=1.0 newton3=enabled\n", mode.c_str(), kCutoff);
  std::printf("# bytesPerParticle=%zu (sizeof(MoleculeLJ)=%zu + %zu SoA arrays) poolMiB=%.1f seed=%zu\n",
              kBytesPerParticle, sizeof(MoleculeType), kSoAArrays, gPoolMiB, gSeed);
  std::printf("# repetitions=%zu iterations=%zu winMargin=%.2f requiredWinRatio=%.2f\n", kRepetitions, kIterations,
              kSortedWinMargin, kRequiredWinRatio);

  std::vector<CellPair> pool;

  if (mode == "itersweep") {
    // Holds the total timed calls fixed and varies how they split between calls-per-region and
    // repetitions, to test whether shortening the inner loop recovers the cold result.
    std::printf("direction,mode,iterations,repetitions,N,ns_per_call_unsorted,ns_per_call_sorted,ratio,win_fraction\n");
    const std::vector<std::pair<size_t, size_t>> splits{{25, 100}, {5, 500}, {1, 2500}};
    for (const auto &direction : directions) {
      const Geometry g = makeGeometry(direction);
      for (const size_t n : {size_t{100}, size_t{160}, size_t{200}}) {
        const size_t seed = gSeed + n;
        for (const auto &[iters, reps] : splits) {
          const RunResult r =
              measure(functor, defaultParticle, g, n, true, /*hot=*/true, pool, targetBytes, seed, iters, reps);
          std::printf("%s,hot,%zu,%zu,%zu,%.1f,%.1f,%.4f,%.3f\n", direction.c_str(), iters, reps, n,
                      r.nsPerCallUnsorted, r.nsPerCallSorted, r.nsPerCallSorted / r.nsPerCallUnsorted, r.winFraction);
          std::fflush(stdout);
        }
        const RunResult c = measure(functor, defaultParticle, g, n, true, /*hot=*/false, pool, targetBytes, seed);
        std::printf("%s,cold,%zu,%zu,%zu,%.1f,%.1f,%.4f,%.3f\n", direction.c_str(), kIterations, kRepetitions, n,
                    c.nsPerCallUnsorted, c.nsPerCallSorted, c.nsPerCallSorted / c.nsPerCallUnsorted, c.winFraction);
        std::fflush(stdout);
      }
    }
    functor.endTraversal(true);
    return 0;
  }

  if (mode == "poolsweep") {
    // Traces the ratio against pool footprint. Footprints are absolute and identical on every
    // machine, so where the sorted path stops losing is a direct read-out of that machine's cache.
    std::printf(
        "direction,N,pool_mib_target,pool_mib_actual,pool_pairs,ns_per_call_unsorted,ns_per_call_sorted,"
        "ratio,win_fraction\n");
    const std::vector<double> targetsMiB{0.0, 1.0, 4.0, 16.0, 64.0, 256.0, 1024.0};
    for (const auto &direction : directions) {
      const Geometry g = makeGeometry(direction);
      for (const size_t n : {size_t{100}, size_t{160}, size_t{200}}) {
        const size_t seed = gSeed + n;
        for (const double t : targetsMiB) {
          const bool hot = t == 0.0;
          const auto bytes = static_cast<size_t>(t * 1024.0 * 1024.0);
          const RunResult r = measure(functor, defaultParticle, g, n, true, hot, pool, bytes, seed);
          std::printf("%s,%zu,%.1f,%.1f,%zu,%.1f,%.1f,%.4f,%.3f\n", direction.c_str(), n, t, r.poolMiB, r.poolPairs,
                      r.nsPerCallUnsorted, r.nsPerCallSorted, r.nsPerCallSorted / r.nsPerCallUnsorted, r.winFraction);
          std::fflush(stdout);
        }
      }
    }
    functor.endTraversal(true);
    return 0;
  }

  if (mode != "full") {
    usage();
    return 1;
  }

  std::printf(
      "direction,mode,N,pool_mib_actual,pool_pairs,ns_per_call_unsorted,ns_per_call_sorted,ratio,"
      "win_fraction\n");
  for (const auto &direction : directions) {
    const Geometry g = makeGeometry(direction);
    for (const bool hot : {true, false}) {
      size_t crossover = 0;
      for (size_t n = 20; n <= 220; n += 10) {
        const RunResult r = measure(functor, defaultParticle, g, n, /*newton3=*/true, hot, pool, targetBytes, gSeed + n);
        std::printf("%s,%s,%zu,%.1f,%zu,%.1f,%.1f,%.4f,%.3f\n", direction.c_str(), hot ? "hot" : "cold", n, r.poolMiB,
                    r.poolPairs, r.nsPerCallUnsorted, r.nsPerCallSorted, r.nsPerCallSorted / r.nsPerCallUnsorted,
                    r.winFraction);
        std::fflush(stdout);
        if (crossover == 0 and r.winFraction >= kRequiredWinRatio) {
          crossover = n;
        }
      }
      std::printf("# CROSSOVER direction=%s mode=%s N=%zu\n", direction.c_str(), hot ? "hot" : "cold", crossover);
      std::fflush(stdout);
    }
  }

  functor.endTraversal(true);
  return 0;
}
