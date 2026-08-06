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
 * Two modes, identical in every other respect (same geometry, same scatter, same alternation,
 * same particle counts, same number of timed calls):
 *   hot  reproduces the production benchmark, kIterations calls on one cell pair.
 *   cold draws each timed call from a pre-built pool of distinct cell pairs large enough that the
 *        pool footprint exceeds last-level cache, so no call reuses the previous call's data.
 *
 * Reports, per direction and per combined particle count, the mean sorted/unsorted time ratio and
 * the fraction of repetitions in which sorting wins by the production margin, for both modes.
 * The crossover under the production decision rule is printed per mode at the end.
 *
 * Not part of the thesis-cited benchmark dataset. Debug branch only.
 */

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <utility>
#include <string>
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

/**
 * Number of distinct cell pairs held in the cold-mode pool. Sized so the pool footprint clears
 * last-level cache. At N=200 combined a pair costs roughly 200 particles * ~200 B ~= 40 kB, so
 * 1000 pairs is ~40 MB, above the 16-32 MB L3 of the Zen 4 parts this runs on. Overridable so the
 * pool can be scaled to a machine with a larger LLC. Must be a multiple of 2 * kIterations so the
 * two halves stay aligned to the repeating split pattern, main() rounds up to enforce that.
 */
size_t gPoolPairs = 1000;

/// Selects the measurement programme. "full" sweeps N in both modes, "itersweep" varies the loop split.
std::string gMode = "full";

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
 * Result of one (direction, N, mode) measurement point.
 */
struct RunResult {
  double meanSortedNs;
  double meanUnsortedNs;
  size_t sortedWins;
};

/**
 * Measures the sorted and unsorted paths at one combined particle count.
 *
 * hot == true reproduces SortingThresholdBenchmark::executeRun() exactly, including the
 * alternation of which path runs first. hot == false performs the same number of timed calls but
 * walks a pool of distinct cell pairs so that no call reuses the previous call's data.
 */
RunResult measure(BenchFunctor &functor, const MoleculeType &defaultParticle, const Geometry &g, size_t numParticles,
                  bool newton3, bool hot, std::vector<CellPair> &pool, size_t iterations = kIterations,
                  size_t repetitions = kRepetitions) {
  BenchCF cellFunctor{functor, kCutoff, autopas::DataLayoutOption::soa, newton3};
  const autopas::SortingThresholdInfoSingle zeroThreshold(0);
  cellFunctor.setSoASortingThresholds(zeroThreshold);
  cellFunctor.setAoSSortingThresholds(zeroThreshold);

  autopas::utils::Timer sortedTimer, unsortedTimer;
  size_t sortedWins = 0;
  std::random_device rd;
  std::mt19937 gen(rd());

  // Cold mode fills a pool whose split pattern repeats with period kIterations, so that the two
  // paths can be handed identical workloads while reading disjoint memory. The sorted path starts
  // half a pool away from the unsorted one, and both walk forward one entry per call, so within a
  // timed region each call touches a pair it has not touched before, exactly as a traversal does.
  const size_t half = pool.size() / 2;
  size_t repOffset = 0;

  if (not hot) {
    std::vector<Split> splitTable(iterations);
    for (auto &s : splitTable) {
      s = makeSplit(numParticles, gen);
    }
    for (size_t p = 0; p < pool.size(); ++p) {
      fillPair(pool[p], g, defaultParticle, splitTable[p % iterations], static_cast<unsigned int>(2 * p));
    }
  }

  for (size_t i = 0; i < repetitions; ++i) {
    CellPair hotPair;
    if (hot) {
      fillPair(hotPair, g, defaultParticle, makeSplit(numParticles, gen), static_cast<unsigned int>(2 * i));
    }

    const auto runOne = [&](bool sorted, autopas::utils::Timer &timer) {
      const std::array<double, 3> dir = sorted ? g.sortingDirection : std::array<double, 3>{0., 0., 0.};
      const size_t base = (repOffset + (sorted ? half : 0)) % pool.size();
      timer.start();
      for (size_t j = 0; j < iterations; ++j) {
        CellPair &pair = hot ? hotPair : pool[(base + j) % pool.size()];
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
    // Advance to pairs neither path has touched yet. Both halves of the pool advance together, so
    // the workload stays matched between the paths.
    repOffset = (repOffset + iterations) % half;
  }

  // Normalise per timed call so points with different iteration counts stay comparable.
  const auto calls = static_cast<double>(repetitions * iterations);
  return {static_cast<double>(sortedTimer.getTotalTime()) / calls,
          static_cast<double>(unsortedTimer.getTotalTime()) / calls, sortedWins};
}

}  // namespace

int main(int argc, char **argv) {
  if (argc > 1) {
    gPoolPairs = std::stoul(argv[1]);
  }
  if (argc > 2) {
    gMode = argv[2];
  }
  const size_t alignment = 2 * kIterations;
  gPoolPairs = ((gPoolPairs + alignment - 1) / alignment) * alignment;

  BenchFunctor functor(kCutoff);
  functor.setParticleProperties(24.0, 1.0);
  functor.initTraversal();
  const MoleculeType defaultParticle({0., 0., 0.}, {0., 0., 0.}, 0, 0);

  const std::vector<std::string> directions{"corner", "edge", "face"};
  std::vector<size_t> counts;
  for (size_t n = 20; n <= 220; n += 10) {
    counts.push_back(n);
  }

  std::printf("# cell edge = cutoff = %.3f, pool pairs = %zu, reps = %zu, iterations/rep = %zu\n", kCutoff, gPoolPairs,
              kRepetitions, kIterations);
  std::printf("# newton3 = enabled\n");
  std::vector<CellPair> pool(gPoolPairs);

  // Iteration-count sweep. Holds the total number of timed calls fixed at 2500 and varies how they
  // are split between calls-per-timed-region and repetitions, to test whether shortening the inner
  // loop recovers the cold result. Reported per timed call so the points are directly comparable.
  if (gMode == "itersweep") {
    std::printf("direction,mode,iterations,repetitions,N,ns_per_call_unsorted,ns_per_call_sorted,ratio,sorted_wins\n");
    const std::vector<std::pair<size_t, size_t>> splits{{25, 100}, {5, 500}, {1, 2500}};
    for (const auto &direction : directions) {
      const Geometry g = makeGeometry(direction);
      for (const size_t n : {size_t{100}, size_t{160}, size_t{200}}) {
        for (const auto &[iters, reps] : splits) {
          const RunResult r = measure(functor, defaultParticle, g, n, true, /*hot=*/true, pool, iters, reps);
          std::printf("%s,hot,%zu,%zu,%zu,%.1f,%.1f,%.4f,%zu\n", direction.c_str(), iters, reps, n, r.meanUnsortedNs,
                      r.meanSortedNs, r.meanSortedNs / r.meanUnsortedNs, r.sortedWins);
          std::fflush(stdout);
        }
        const RunResult c = measure(functor, defaultParticle, g, n, true, /*hot=*/false, pool, kIterations,
                                    kRepetitions);
        std::printf("%s,cold,%zu,%zu,%zu,%.1f,%.1f,%.4f,%zu\n", direction.c_str(), kIterations, kRepetitions, n,
                    c.meanUnsortedNs, c.meanSortedNs, c.meanSortedNs / c.meanUnsortedNs, c.sortedWins);
        std::fflush(stdout);
      }
    }
    functor.endTraversal(true);
    return 0;
  }

  std::printf("direction,mode,N,mean_unsorted_ns,mean_sorted_ns,ratio,sorted_wins\n");
  for (const auto &direction : directions) {
    const Geometry g = makeGeometry(direction);
    for (const bool hot : {true, false}) {
      size_t crossover = 0;
      for (const size_t n : counts) {
        const RunResult r = measure(functor, defaultParticle, g, n, /*newton3=*/true, hot, pool);
        const double ratio = r.meanSortedNs / r.meanUnsortedNs;
        std::printf("%s,%s,%zu,%.1f,%.1f,%.4f,%zu\n", direction.c_str(), hot ? "hot" : "cold", n, r.meanUnsortedNs,
                    r.meanSortedNs, ratio, r.sortedWins);
        std::fflush(stdout);
        if (crossover == 0 and
            static_cast<double>(r.sortedWins) / static_cast<double>(kRepetitions) >= kRequiredWinRatio) {
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
