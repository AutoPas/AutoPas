/**
 * @file SortingThresholdBenchmark.h
 * @date 27.06.2026
 * @author hmeyran
 */

#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <random>
#include <string_view>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Determines the per-Newton3-state, per-direction-type particle-count threshold at which using a sorted cell-pair
 * interaction path becomes faster than the unsorted path, separately for the AoS path (SortedCellView vs. the
 * unsorted double loop) and the SoA path (SoAFunctorPairSorted vs. SoAFunctorPair).
 *
 * The AoS half is benchmarked for any pairwise functor. The SoA half is only benchmarked for functors with
 * Functor_T::supportsSoASorting == true (see runBenchmark()).
 *
 * Stores one threshold per Newton3 state and direction type, indexed as [newton3][direction]:
 *   [0]   → Newton3 disabled, [1] → Newton3 enabled
 *   [][0] → Corner (three non-zero components)
 *   [][1] → Edge   (two non-zero components)
 *   [][2] → Face   (one non-zero component)
 *
 * Newton3 on/off is benchmarked separately because it changes how much work the sorted/unsorted paths do per call
 * (symmetric force application vs. one-sided passes), which can shift the break-even point.
 *
 * The search performed in runSearch() is deliberately biased towards a conservative (higher) threshold, to avoid noisy
 * measurements influencing the threshold to be too low, causing slowdowns.
 */
class SortingThresholdBenchmark {
 public:
  /**
   * Return all per-Newton3-state, per-direction-type SoA sorting thresholds.
   * Only meaningful once hasRunSoA() is true.
   * @return Copy of the internal threshold array, indexed as [newton3][direction] (see class documentation).
   */
  std::array<std::array<size_t, 3>, 2> getSoAThresholds() const { return _soaThresholds; }

  /**
   * Return all per-Newton3-state, per-direction-type AoS pair-sorting thresholds.
   * Only meaningful once hasRunAoS() is true.
   * @return Copy of the internal threshold array, indexed as [newton3][direction] (see class documentation).
   */
  std::array<std::array<size_t, 3>, 2> getAoSThresholds() const { return _aosThresholds; }

  /**
   * Returns whether runBenchmark<..., true>() (the SoA half) has already been called.
   * @return True if the SoA benchmark has run.
   */
  [[nodiscard]] bool hasRunSoA() const { return _hasRunSoA; }

  /**
   * Returns whether runBenchmark<..., false>() (the AoS half) has already been called.
   * @return True if the AoS benchmark has run.
   */
  [[nodiscard]] bool hasRunAoS() const { return _hasRunAoS; }

  /**
   * Runs the micro-benchmark for sorting thresholds, sweeping both Newton3 states and all three direction types,
   * and storing the resulting thresholds.
   * Templated on which path to benchmark so callers can trigger each path
   * independently and only instantiate this for the functor/particle combinations they actually need.
   * - UseSoA == true: benchmarks the SoA path and sets hasRunSoA() to true. Callers should only instantiate this
   *   for functors with Functor_T::supportsSoASorting == true.
   * - UseSoA == false: benchmarks the AoS path and sets hasRunAoS() to true. Internally
   *   guarded by whether Particle_T is constructible the way this benchmark needs (position, velocity, id):.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @tparam UseSoA Whether to benchmark the SoA path (true) or the AoS path (false).
   * @param functor Functor instance used to drive the benchmark cells.
   */
  template <class Functor_T, class Particle_T, bool UseSoA>
  void runBenchmark(Functor_T &functor) {
    if constexpr (UseSoA) {
      for (size_t n3 = 0; n3 < 2; n3++) {
        for (size_t layout = 0; layout < 3; layout++) {
          _soaThresholds[n3][layout] = runSearch<Functor_T, Particle_T, true>(functor, layout, static_cast<bool>(n3));
        }
      }
      _hasRunSoA = true;
    } else {
      if constexpr (std::constructible_from<Particle_T, std::array<double, 3>, std::array<double, 3>, size_t>) {
        for (size_t n3 = 0; n3 < 2; n3++) {
          for (size_t layout = 0; layout < 3; layout++) {
            _aosThresholds[n3][layout] =
                runSearch<Functor_T, Particle_T, false>(functor, layout, static_cast<bool>(n3));
          }
        }
      } else {
        AutoPasLog(WARN,
                   "SortingThresholdBenchmark: Particle type is not constructible from (position, velocity, id), "
                   "skipping the AoS sorting-threshold benchmark and using the default thresholds instead.");
      }
      // Here so that this benchmark doesn't get executed every iteration if the particle doesn't support threshold
      // benchmarking
      _hasRunAoS = true;
    }
  }

 private:
  /**
   * Set to true by runBenchmark<..., true>() once the SoA benchmark has completed.
   */
  bool _hasRunSoA{false};

  /**
   * Set to true by runBenchmark<..., false>() once the AoS benchmark has completed.
   */
  bool _hasRunAoS{false};

  /**
   * Per-Newton3-state, per-direction-type SoA sorting threshold values, initialized to the compile-time default.
   * Indexed as [newton3][direction] (see class documentation for the indexing convention).
   */
  std::array<std::array<size_t, 3>, 2> _soaThresholds{{{100, 100, 100}, {100, 100, 100}}};

  /**
   * Per-Newton3-state, per-direction-type AoS pair-sorting threshold values, initialized to the compile-time
   * default. Indexed as [newton3][direction] (see class documentation for the indexing convention).
   */
  std::array<std::array<size_t, 3>, 2> _aosThresholds{{{8, 8, 8}, {8, 8, 8}}};

  /**
   * Human-readable names for the direction-type index (0=Corner, 1=Edge, 2=Face), used only for log messages.
   */
  static constexpr std::array<std::string_view, 3> _layoutNames{"Corner", "Edge", "Face"};

  /**
   * Human-readable names for the Newton3-state index (0=disabled, 1=enabled), used only for log messages.
   */
  static constexpr std::array<std::string_view, 2> _newton3Names{"N3off", "N3on"};

  /**
   * Number of timed calls per repetition to get stable measurement.
   * @todo: Adjust this based on input size.
   */
  const size_t _iterations = 25;

  /**
   * Number of independent measurement repetitions per particle count.
   * @todo: Try out different numbers of repetitions.
   */
  const size_t _repetitions = 25;
  /**
   * Upper bound on the particle count searched by the binary search for the SoA path.
   */
  const size_t _maxSoAParticles = 250;

  /**
   * Upper bound on the particle count searched by the binary search for the AoS path.
   * Smaller than _maxSoAParticles because the AoS sorted path (SortedCellView) is expected to pay off at
   * substantially lower particle counts than the SoA sorted path.
   */
  const size_t _maxAoSParticles = 50;

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
   * Factor of particles that get randomly split between the two cells to not always have clean 50/50 splits.
   */
  const double _scatterFactor = 0.2;

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
   * Measures the mean per-repetition time for the sorted and unsorted pair interaction paths at a given particle
   * count for one direction type and Newton3 state, for either the AoS or the SoA layout.
   *
   * Particles are regenerated each repetition so the sort always operates on fresh random data.
   * The sorted path is controlled by passing the layout-specific sortingDirection; the unsorted
   * path is forced by passing a zero direction (CellFunctor skips sorting when direction is zero).
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @tparam useSoA Selects the benchmarked data layout: true measures the SoA path, false the AoS path.
   * @param functor Functor instance used to drive the benchmark.
   * @param layout Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param numParticles Number of particles placed in each of the two cells.
   * @param useNewton3 Whether the benchmarked CellFunctor should apply Newton3.
   * @return The number of repetitions (out of _repetitions) in which the sorted path was measured as faster than
   * the unsorted path by at least _sortedWinMarginFraction.
   */
  template <class Functor_T, class Particle_T, bool useSoA>
  size_t executeRun(Functor_T &functor, size_t layout, size_t numParticles, bool useNewton3) {
    using BenchCell = FullParticleCell<Particle_T>;
    using BenchCF = internal::CellFunctor<BenchCell, Functor_T, false>;

    const Particle_T defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
    const double cutoff = functor.getCutoff();
    const double invSqrt3 = 1. / sqrt(3.);
    const double invSqrt2 = 1. / sqrt(2.);

    size_t numParticlesEqDistr = static_cast<size_t>(std::ceil(numParticles * (1. - _scatterFactor)));

    size_t numParticlesScatter = numParticles - numParticlesEqDistr;

    size_t numParticlesCell1 = numParticlesEqDistr / 2;
    size_t numParticlesCell2 = numParticlesEqDistr / 2 + numParticlesEqDistr % 2;

    // Prepare for random scattering.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distrib(0, numParticlesScatter);

    // Setup Cell Layout
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
      default:
        utils::ExceptionHandler::exception("Layout {} is not a valid/supported layout!", layout);
    }

    utils::Timer sortedTimer, unsortedTimer;
    BenchCF cellFunctor{functor, cutoff, useSoA ? DataLayoutOption::soa : DataLayoutOption::aos, useNewton3};
    size_t sortedWins = 0;
    // Set to 0 so whether sorting happens is controlled entirely through the sorting direction.
    cellFunctor.setSoASortingThreshold(0);
    cellFunctor.setAoSSortingThreshold(0);

    // For SoA Forces won't accumulate, as SoAExtractor is never called, for AoS they will as AoS modifies the values
    // directly.
    auto measureUnsorted = [&]() {
      unsortedTimer.start();
      for (size_t j = 0; j < _iterations; j++) {
        if constexpr (useSoA) {
          functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
          functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
        }
        // A sorting direction of (0, 0, 0) disables sorting.
        cellFunctor.processCellPair(cell1, cell2, {0., 0., 0.});
      }
      return unsortedTimer.stop();
    };
    auto measureSorted = [&]() {
      sortedTimer.start();
      for (size_t j = 0; j < _iterations; j++) {
        if constexpr (useSoA) {
          functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
          functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
        }
        cellFunctor.processCellPair(cell1, cell2, sortingDirection);
      }
      return sortedTimer.stop();
    };

    for (size_t i = 0; i < _repetitions; i++) {
      cell1.clear();
      cell2.clear();

      // Generate the number of scattered particles for cell1 and cell2
      size_t toAddCell1 = distrib(gen);
      // Vary the seed per repetition per cell so each repetition samples a fresh particle
      // layout instead of repeatedly timing the exact same configuration.
      fillWithRandomParticles(cell1, defaultParticle, cell1Low, cell1High, numParticlesCell1 + toAddCell1,
                              static_cast<unsigned int>(2 * i));
      fillWithRandomParticles(cell2, defaultParticle, cell2Low, cell2High,
                              numParticlesCell2 + numParticlesScatter - toAddCell1,
                              static_cast<unsigned int>(2 * i + 1));

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

      AutoPasLog(TRACE, "SortingThresholdBenchmark rep {}/{} layout={} n={}: unsorted={}ns sorted={}ns", i + 1,
                 _repetitions, _layoutNames[layout], numParticles, unsortedDelta, sortedDelta);

      // A repetition only counts as a "sorted win" if it clears the margin: see _sortedWinMarginFraction.
      if (static_cast<double>(sortedDelta) < static_cast<double>(unsortedDelta) * (1. - _sortedWinMarginFraction)) {
        ++sortedWins;
      }
    }

    const long meanSorted = sortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    const long meanUnsorted = unsortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    AutoPasLog(TRACE, "SortingThresholdBenchmark layout={} n={}: mean unsorted={}ns mean sorted={}ns sortedWins={}/{}",
               _layoutNames[layout], numParticles, meanUnsorted, meanSorted, sortedWins, _repetitions);
    return sortedWins;
  }

  /**
   * Binary-searches over particle count to find the smallest n at which the sorted path is faster
   * than the unsorted path for a given direction type, for either the AoS or the SoA layout.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @tparam useSoA Selects the benchmarked data layout: true = SoA, false = AoS.
   * @param functor Functor instance used to drive the benchmark.
   * @param layout Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param useNewton3 Whether the benchmarked CellFunctor should apply Newton3.
   * @return Smallest particle count at which sorted beats unsorted, or the upper search bound if never.
   */
  template <class Functor_T, class Particle_T, bool useSoA>
  size_t runSearch(Functor_T &functor, size_t layout, bool useNewton3) {
    size_t lowCount = 0;
    // Use a smaller max particle count for AoS, since the AoS sorted path is expected to pay off at a lower
    // particle count than the SoA sorted path (see _maxAoSParticles).
    size_t highCount = useSoA ? _maxSoAParticles : _maxAoSParticles;

    while (lowCount < highCount) {
      size_t mid = lowCount + (highCount - lowCount) / 2;

      const auto outcome = executeRun<Functor_T, Particle_T, useSoA>(functor, layout, mid, useNewton3);
      const double winRatio = static_cast<double>(outcome) / static_cast<double>(_repetitions);

      // Conservative decision rule: only accept "sorted wins" once a clear majority of repetitions
      // agree by a clear margin (see _sortedWinMarginFraction and _requiredSortedWinRatio).
      if (winRatio >= _requiredSortedWinRatio) {
        highCount = mid;
        AutoPasLog(TRACE, "SortingThresholdBenchmark search {} layout={} n={}: sorted won {}/{} reps → high={}",
                   _newton3Names[useNewton3], _layoutNames[layout], mid, outcome, _repetitions, highCount);
      } else {
        lowCount = mid + 1;
        AutoPasLog(TRACE, "SortingThresholdBenchmark search {} layout={} n={}: sorted won only {}/{} reps → low={}",
                   _newton3Names[useNewton3], _layoutNames[layout], mid, outcome, _repetitions, lowCount);
      }
    }
    AutoPasLog(DEBUG, "SortingThresholdBenchmark {} layout={} threshold={}", _newton3Names[useNewton3],
               _layoutNames[layout], lowCount);
    return lowCount;
  }
};

}  // namespace autopas
