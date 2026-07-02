/**
 * @file SortingThresholdBenchmark.h
 * @date 27.06.2026
 * @author hmeyran
 */

#pragma once

#include <array>
#include <cmath>
#include <random>
#include <string_view>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Determines the per-direction-type particle-count threshold at which using the sorted SoA path
 * (SoAFunctorPairSorted) becomes faster than the unsorted path (SoAFunctorPair).
 *
 * Stores one threshold per direction type, indexed by the number of zero components in the
 * sortingDirection vector:
 *   0 → Corner (three non-zero components)
 *   1 → Edge   (two non-zero components)
 *   2 → Face   (one non-zero component)
 *
 * The search performed in runSearch() is deliberately biased towards a conservative (higher) threshold, to avoid noisy
 * measurements to influence the threshold to low, causing slow downs.
 */
class SortingThresholdBenchmark {
 public:
  /**
   * Return all three per-direction-type thresholds.
   * Index 0 = Corner (zero zeros), 1 = Edge (one zero), 2 = Face (two zeros).
   * @return Copy of the internal threshold array.
   */
  std::array<size_t, 3> getThresholds() const { return _thresholds; }

  /**
   * Returns whether runBenchmark() has already been called.
   * @return True if the benchmark has run.
   */
  [[nodiscard]] bool hasRun() const { return _hasRun; }

  /**
   * Runs the micro-benchmark for all three direction types and stores the resulting thresholds.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @param functor Functor instance used to drive the benchmark cells.
   */
  template <class Functor_T, class Particle_T>
  void runBenchmark(Functor_T &functor) {
    for (size_t layout = 0; layout < 3; layout++) {
      _thresholds[layout] = runSearch<Functor_T, Particle_T>(functor, layout);
    }
    _hasRun = true;
  }

 private:
  /**
   * Set to true by runBenchmark() once the benchmark has completed.
   */
  bool _hasRun{false};

  /**
   * Per-direction-type threshold values, initialized to the compile-time default.
   * Index = number of zero components in sortingDirection (0=Corner, 1=Edge, 2=Face).
   */
  std::array<size_t, 3> _thresholds{25, 25, 25};

  /**
   * Human-readable names for the direction-type index (0=Corner, 1=Edge, 2=Face), used only for log messages.
   */
  static constexpr std::array<std::string_view, 3> _layoutNames{"Corner", "Edge", "Face"};

  /**
   * Number of timed calls per repetition; amortizes timer overhead within a single rep.
   * @todo: Adjust this based on input size.
   */
  const size_t _iterations = 100;

  /**
   * Number of independent measurement repetitions per particle count; the mean is taken over these.
   * @todo: Try out different numbers of repetitions.
   */
  const size_t _repetitions = 25;
  /**
   * Upper bound on the particle count searched by the binary search.
   */
  const size_t _maxParticles = 500;

  /**
   * Minimum relative speed-up the sorted path must show over the unsorted path within a single repetition
   * for that repetition to count as a clear win for sorting.
   * @todo: Adjust this based on how noisy the target hardware is.
   */
  const double _sortedWinMarginFraction = 0.05;

  /**
   * Minimum fraction of repetitions within a single executeRun() call that must count as a clear win for
   * sorting (see _sortedWinMarginFraction) before runSearch() accepts "sorted is faster" for the tested
   * particle count.
   * @todo: Adjust this based on how noisy the target hardware is.
   */
  const double _requiredSortedWinRatio = 0.7;

  /**
   * Fills a cell with numParticles copies of defaultParticle at independently uniform-random positions within
   * [boxLow, boxHigh]. A minimal stand-in for autopasTools::generators::UniformGenerator::fillWithParticles() so
   * that this core-library header does not need to depend on the autopasTools target.
   * @tparam Cell_T Cell type; must support addParticle().
   * @tparam Particle_T Particle type.
   * @param cell Cell to fill.
   * @param defaultParticle Template particle whose non-positional properties are copied.
   * @param boxLow Lower corner of the sampling box.
   * @param boxHigh Upper corner of the sampling box.
   * @param numParticles Number of particles to generate.
   * @param seed Seed for the random engine.
   */
  template <class Cell_T, class Particle_T>
  static void fillWithRandomParticles(Cell_T &cell, const Particle_T &defaultParticle,
                                      const std::array<double, 3> &boxLow, const std::array<double, 3> &boxHigh,
                                      size_t numParticles, unsigned int seed = 42) {
    std::mt19937 generator(seed);
    std::array<std::uniform_real_distribution<double>, 3> dist{
        std::uniform_real_distribution<double>(boxLow[0], boxHigh[0]),
        std::uniform_real_distribution<double>(boxLow[1], boxHigh[1]),
        std::uniform_real_distribution<double>(boxLow[2], boxHigh[2])};
    for (size_t i = 0; i < numParticles; ++i) {
      Particle_T particle(defaultParticle);
      particle.setR({dist[0](generator), dist[1](generator), dist[2](generator)});
      particle.setID(i);
      particle.setOwnershipState(autopas::OwnershipState::owned);
      cell.addParticle(particle);
    }
  }

  /**
   * Measures the mean per-repetition time for the sorted and unsorted SoA pair interaction paths
   * at a given particle count for one direction type.
   *
   * Particles are regenerated each repetition so the sort always operates on fresh random data.
   * The sorted path is controlled by passing the layout-specific sortingDirection; the unsorted
   * path is forced by passing a zero direction (CellFunctor skips sorting when direction is zero).
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @param functor Functor instance used to drive the benchmark.
   * @param layout Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param numParticles Number of particles placed in each of the two cells.
   * @return Pooled mean sorted/unsorted runtimes and the per-repetition sorted-win count.
   */
  template <class Functor_T, class Particle_T>
  size_t executeRun(Functor_T &functor, size_t layout, size_t numParticles) {
    using BenchCell = FullParticleCell<Particle_T>;
    using BenchCF = internal::CellFunctor<BenchCell, Functor_T, false>;

    const Particle_T defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
    const double cutoff = functor.getCutoff();
    const double invSqrt3 = 1. / sqrt(3.);
    const double invSqrt2 = 1. / sqrt(2.);
    BenchCF cellFunctor{functor, cutoff, DataLayoutOption::soa, false};
    // Set to 0 so whether sorting happens is controlled entirely through the sorting direction.
    cellFunctor.setSoASortingThreshold(0);
    BenchCell cell1, cell2;

    std::array cell1Low = {0., 0., 0.};
    std::array cell1High = {cutoff, cutoff, cutoff};

    std::array cell2Low = {0., 0., 0.};
    std::array cell2High = {cutoff, cutoff, cutoff};

    std::array sortingDirection = {0., 0., 0.};

    switch (layout) {
      case 0:
        cell2Low = {cutoff, cutoff, cutoff};
        cell2High = {2. * cutoff, 2. * cutoff, 2. * cutoff};
        sortingDirection = {invSqrt3, invSqrt3, invSqrt3};
        break;
      case 1:
        cell2Low = {cutoff, cutoff, 0.};
        cell2High = {2. * cutoff, 2. * cutoff, cutoff};
        sortingDirection = {invSqrt2, invSqrt2, 0.};
        break;
      case 2:
        cell2Low = {cutoff, 0., 0.};
        cell2High = {2 * cutoff, cutoff, cutoff};
        sortingDirection = {1, 0., 0.};
        break;
    }
    utils::Timer sortedTimer, unsortedTimer;
    size_t sortedWins = 0;
    for (size_t i = 0; i < _repetitions; i++) {
      cell1.clear();
      cell2.clear();

      // Vary the seed per repetition per cell so each repetition samples a fresh particle
      // layout instead of repeatedly timing the exact same configuration.
      fillWithRandomParticles(cell1, defaultParticle, cell1Low, cell1High, numParticles,
                              static_cast<unsigned int>(2 * i));
      fillWithRandomParticles(cell2, defaultParticle, cell2Low, cell2High, numParticles,
                              static_cast<unsigned int>(2 * i + 1));
      functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
      functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);

      auto measureUnsorted = [&]() {
        unsortedTimer.start();
        for (size_t j = 0; j < _iterations; j++) {
          // A sorting direction of (0, 0, 0) disables sorting.
          cellFunctor.processCellPair(cell1, cell2, {0., 0., 0.});
        }
        return unsortedTimer.stop();
      };
      auto measureSorted = [&]() {
        sortedTimer.start();
        for (size_t j = 0; j < _iterations; j++) {
          cellFunctor.processCellPair(cell1, cell2, sortingDirection);
        }
        return sortedTimer.stop();
      };

      long unsortedDelta = 0;
      long sortedDelta = 0;
      // Alternate which path is measured first. Whichever path runs second inherits warm caches and a settled
      // branch predictor from the first, which would otherwise make it look systematically faster than it is.
      if (i % 2 == 0) {
        unsortedDelta = measureUnsorted();
        sortedDelta = measureSorted();
      } else {
        sortedDelta = measureSorted();
        unsortedDelta = measureUnsorted();
      }

      AutoPasLog(DEBUG, "SortingThresholdBenchmark rep {}/{} layout={} n={}: unsorted={}ns sorted={}ns", i + 1,
                 _repetitions, _layoutNames[layout], numParticles, unsortedDelta, sortedDelta);

      // A repetition only counts as a "sorted win" if it clears the margin: see _sortedWinMarginFraction.
      if (static_cast<double>(sortedDelta) < static_cast<double>(unsortedDelta) * (1. - _sortedWinMarginFraction)) {
        ++sortedWins;
      }
    }

    const long meanSorted = sortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    const long meanUnsorted = unsortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    AutoPasLog(INFO, "SortingThresholdBenchmark layout={} n={}: mean unsorted={}ns mean sorted={}ns sortedWins={}/{}",
               _layoutNames[layout], numParticles, meanUnsorted, meanSorted, sortedWins, _repetitions);
    return sortedWins;
  }

  /**
   * Binary-searches over particle count to find the smallest n at which the sorted path is faster
   * than the unsorted path for a given direction type.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @param functor Functor instance used to drive the benchmark.
   * @param layout Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @return Smallest particle count at which sorted beats unsorted, or the upper search bound if never.
   */
  template <class Functor_T, class Particle_T>
  size_t runSearch(Functor_T &functor, size_t layout) {
    size_t lowCount = 0;
    size_t highCount = _maxParticles;

    while (lowCount < highCount) {
      size_t mid = lowCount + (highCount - lowCount) / 2;

      const auto outcome = executeRun<Functor_T, Particle_T>(functor, layout, mid);
      const double winRatio = static_cast<double>(outcome) / static_cast<double>(_repetitions);

      // Conservative decision rule: only accept "sorted wins" once a clear majority of repetitions
      // agree by a clear margin (see _sortedWinMarginFraction and _requiredSortedWinRatio).
      if (winRatio >= _requiredSortedWinRatio) {
        highCount = mid;
        AutoPasLog(DEBUG, "SortingThresholdBenchmark search layout={} n={}: sorted won {}/{} reps → high={}",
                   _layoutNames[layout], mid, outcome, _repetitions, highCount);
      } else {
        lowCount = mid + 1;
        AutoPasLog(DEBUG, "SortingThresholdBenchmark search layout={} n={}: sorted won only {}/{} reps → low={}",
                   _layoutNames[layout], mid, outcome, _repetitions, lowCount);
      }
    }
    AutoPasLog(INFO, "SortingThresholdBenchmark layout={} threshold={}", _layoutNames[layout], lowCount);
    return lowCount;
  }
};

}  // namespace autopas
