/**
 * @file SortingThresholdBenchmark.h
 * @date 27.06.2026
 * @author hmeyran
 */

#pragma once

#include <array>
#include <cmath>
#include <memory>
#include <random>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "autopas/utils/SortingThresholdInfo2B.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {

/**
 * Determines the per-Newton3-state, per-direction-type particle-count threshold at which using a sorted cell-pair
 * interaction path becomes faster than the unsorted path, separately for the AoS path (SortedCellView vs. the
 * unsorted double loop) and the SoA path (SoAFunctorPairSorted vs. SoAFunctorPair).
 *
 * The AoS half is benchmarked for any pairwise functor. The SoA half is only benchmarked for functors with
 * Functor_T::supportsSoASorting == true (see runSoABenchmark()).
 *
 * Stores one threshold per Newton3 state and direction type.
 *
 * Newton3 on/off is benchmarked separately because it changes how much work the sorted/unsorted paths do per call
 * (symmetric force application vs. one-sided passes), which can shift the break-even point.
 *
 * The search performed in runSearch() is deliberately biased towards a conservative (higher) threshold, to avoid noisy
 * measurements influencing the threshold to be too low. A too low threshold would actively regress performance while a
 * too high threshold would at worst still perform the same as baseline. One goal of these thresholds is to avoid any
 * performance regression.
 */
class SortingThresholdBenchmark {
 public:
  /**
   * Return all per-Newton3-state, per-direction-type SoA sorting thresholds.
   * Only meaningful once hasRunSoA() is true.
   * The returned pointer is stable (the same instance is reused every call) once runSoABenchmark() has completed,
   * so callers can cheaply re-apply it every iteration (e.g. to a container) without reallocating.
   * @return Shared pointer to the internal threshold values, indexed as [newton3][direction] (see class
   * documentation).
   */
  std::shared_ptr<const SortingThresholdInfo2B> getSoAThreshold() const { return _soaThresholds; }

  /**
   * Return all per-Newton3-state, per-direction-type AoS pair-sorting thresholds.
   * Only meaningful once hasRunAoS() is true.
   * The returned pointer is stable (the same instance is reused every call) once runAoSBenchmark() has completed,
   * so callers can cheaply re-apply it every iteration (e.g. to a container) without reallocating.
   * @return Shared pointer to the internal threshold values, indexed as [newton3][direction] (see class
   * documentation).
   */
  std::shared_ptr<const SortingThresholdInfo2B> getAoSThreshold() const { return _aosThresholds; }

  /**
   * Returns whether runSoABenchmark() has already been called.
   * @return True if the SoA benchmark has run.
   */
  [[nodiscard]] bool hasRunSoA() const { return _hasRunSoA; }

  /**
   * Returns whether runAoSBenchmark() has already been called.
   * @return True if the AoS benchmark has run.
   */
  [[nodiscard]] bool hasRunAoS() const { return _hasRunAoS; }

  /**
   * Runs the micro-benchmark for the SoA sorting threshold, sweeping both Newton3 states and all three direction
   * types, and storing the resulting thresholds. Sets hasRunSoA() to true.
   * Callers should only instantiate this for functors with Functor_T::supportsSoASorting == true.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @param functor Functor instance used to drive the benchmark cells.
   * @param defaultParticle Template particle (e.g. sampled from the live simulation) whose non-positional
   * properties are copied onto every particle generated for the benchmark cells.
   */
  template <class Functor_T, class Particle_T>
  void runSoABenchmark(Functor_T &functor, const Particle_T &defaultParticle) {
    std::array<std::array<size_t, 3>, 2> thresholds;
    for (const auto &n3 : Newton3Option::getAllOptions()) {
      for (const auto &layout : CellLayoutOption::getAllOptions()) {
        thresholds[n3][layout] = runSearch<Functor_T, Particle_T, true>(functor, defaultParticle, layout, n3);
      }
    }
    // thresholds[n3][layout] is indexed by CellLayoutOption (corner=0, edge=1, face=2), while the constructor
    // below takes named (Face, Edge, Corner) parameters per Newton3 state, so the layout index is reversed here.
    _soaThresholds = std::make_shared<const SortingThresholdInfo2B>(
        thresholds[0][2], thresholds[0][1], thresholds[0][0], thresholds[1][2], thresholds[1][1], thresholds[1][0]);
    _hasRunSoA = true;
  }

  /**
   * Runs the micro-benchmark for the AoS sorting threshold, sweeping both Newton3 states and all three direction
   * types, and storing the resulting thresholds. Sets hasRunAoS() to true.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @param functor Functor instance used to drive the benchmark cells.
   * @param defaultParticle Template particle (e.g. sampled from the live simulation) whose non-positional
   * properties are copied onto every particle generated for the benchmark cells.
   */
  template <class Functor_T, class Particle_T>
  void runAoSBenchmark(Functor_T &functor, const Particle_T &defaultParticle) {
    std::array<std::array<size_t, 3>, 2> thresholds;
    for (const auto &n3 : Newton3Option::getAllOptions()) {
      for (const auto &layout : CellLayoutOption::getAllOptions()) {
        thresholds[n3][layout] = runSearch<Functor_T, Particle_T, false>(functor, defaultParticle, layout, n3);
      }
    }
    // thresholds[n3][layout] is indexed by CellLayoutOption (corner=0, edge=1, face=2), while the constructor
    // below takes named (Face, Edge, Corner) parameters per Newton3 state, so the layout index is reversed here.
    _aosThresholds = std::make_shared<const SortingThresholdInfo2B>(
        thresholds[0][2], thresholds[0][1], thresholds[0][0], thresholds[1][2], thresholds[1][1], thresholds[1][0]);
    _hasRunAoS = true;
  }

 private:
  /**
   * Set to true by runSoABenchmark() once the SoA benchmark has completed.
   */
  bool _hasRunSoA{false};

  /**
   * Set to true by runAoSBenchmark() once the AoS benchmark has completed.
   */
  bool _hasRunAoS{false};

  /**
   * Per-Newton3-state, per-direction-type SoA sorting threshold values, initialized to the compile-time default.
   * Held as a shared_ptr so that getSoAThreshold() can hand out the same instance on every call once
   * runSoABenchmark() has completed, instead of forcing callers to reallocate on every access.
   */
  std::shared_ptr<const SortingThresholdInfo2B> _soaThresholds{std::make_shared<const SortingThresholdInfo2B>(100)};

  /**
   * Per-Newton3-state, per-direction-type AoS pair-sorting threshold values, initialized to the compile-time
   * default.
   * Held as a shared_ptr so that getAoSThreshold() can hand out the same instance on every call once
   * runAoSBenchmark() has completed, instead of forcing callers to reallocate on every access.
   */
  std::shared_ptr<const SortingThresholdInfo2B> _aosThresholds{std::make_shared<const SortingThresholdInfo2B>(8)};

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
   * @param defaultParticle Template particle whose non-positional properties are copied onto every particle
   * generated for the benchmark cells.
   * @param layout Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param numParticles Number of particles placed in each of the two cells.
   * @param newton3 Whether the benchmarked CellFunctor should apply Newton3.
   * @return The number of repetitions (out of _repetitions) in which the sorted path was measured as faster than
   * the unsorted path by at least _sortedWinMarginFraction.
   */
  template <class Functor_T, class Particle_T, bool useSoA>
  size_t executeRun(Functor_T &functor, const Particle_T &defaultParticle, CellLayoutOption layout, size_t numParticles,
                    Newton3Option newton3) {
    using BenchCell = FullParticleCell<Particle_T>;
    using BenchCF = internal::CellFunctor<BenchCell, Functor_T, false>;

    const double cutoff = functor.getCutoff();
    const double invSqrt3 = 1. / sqrt(3.);
    const double invSqrt2 = 1. / sqrt(2.);

    auto numParticlesEqDistr = static_cast<size_t>(std::ceil(numParticles * (1. - _scatterFactor)));
    size_t numParticlesScatter = numParticles - numParticlesEqDistr + numParticlesEqDistr % 2;
    size_t numParticlesPerCell = numParticlesEqDistr / 2;

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
      case CellLayoutOption::corner:
        cell2Low = {cutoff, cutoff, cutoff};
        cell2High = {2. * cutoff, 2. * cutoff, 2. * cutoff};
        sortingDirection = {invSqrt3, invSqrt3, invSqrt3};
        break;
      case CellLayoutOption::edge:
        cell2Low = {cutoff, cutoff, 0.};
        cell2High = {2. * cutoff, 2. * cutoff, cutoff};
        sortingDirection = {invSqrt2, invSqrt2, 0.};
        break;
      case CellLayoutOption::face:
        cell2Low = {cutoff, 0., 0.};
        cell2High = {2 * cutoff, cutoff, cutoff};
        sortingDirection = {1, 0., 0.};
        break;
      default:
        utils::ExceptionHandler::exception("Layout {} is not a valid/supported layout!", layout);
    }

    utils::Timer sortedTimer, unsortedTimer;
    BenchCF cellFunctor{functor, cutoff, useSoA ? DataLayoutOption::soa : DataLayoutOption::aos,
                        newton3 == Newton3Option::enabled};
    size_t sortedWins = 0;
    // Set to 0 so whether sorting happens is controlled entirely through the sorting direction.
    const SortingThresholdInfoSingle zeroThreshold(0);
    cellFunctor.setSoASortingThresholds(zeroThreshold);
    cellFunctor.setAoSSortingThresholds(zeroThreshold);

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
      generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, cell1Low, cell1High,
                                                      numParticlesPerCell + toAddCell1,
                                                      static_cast<unsigned int>(2 * i));
      generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, cell2Low, cell2High,
                                                      numParticlesPerCell + numParticlesScatter - toAddCell1,
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
                 _repetitions, layout, numParticles, unsortedDelta, sortedDelta);

      // A repetition only counts as a "sorted win" if it clears the margin: see _sortedWinMarginFraction.
      if (static_cast<double>(sortedDelta) < static_cast<double>(unsortedDelta) * (1. - _sortedWinMarginFraction)) {
        ++sortedWins;
      }
    }

    const long meanSorted = sortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    const long meanUnsorted = unsortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    AutoPasLog(TRACE, "SortingThresholdBenchmark layout={} n={}: mean unsorted={}ns mean sorted={}ns sortedWins={}/{}",
               layout, numParticles, meanUnsorted, meanSorted, sortedWins, _repetitions);
    return sortedWins;
  }

  /**
   * Binary-searches over particle count to find the smallest n at which the sorted path is faster
   * than the unsorted path for a given direction type, for either the AoS or the SoA layout.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @tparam useSoA Selects the benchmarked data layout: true = SoA, false = AoS.
   * @param functor Functor instance used to drive the benchmark.
   * @param defaultParticle Template particle whose non-positional properties are copied onto every particle
   * generated for the benchmark cells.
   * @param layout Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param newton3 Whether the benchmarked CellFunctor should apply Newton3.
   * @return Smallest particle count at which sorted beats unsorted, or the upper search bound if never.
   */
  template <class Functor_T, class Particle_T, bool useSoA>
  size_t runSearch(Functor_T &functor, const Particle_T &defaultParticle, CellLayoutOption layout,
                   Newton3Option newton3) {
    size_t lowCount = 0;
    size_t highCount = useSoA ? _maxSoAParticles : _maxAoSParticles;

    while (lowCount < highCount) {
      size_t mid = lowCount + (highCount - lowCount) / 2;

      const auto outcome = executeRun<Functor_T, Particle_T, useSoA>(functor, defaultParticle, layout, mid, newton3);
      const double winRatio = static_cast<double>(outcome) / static_cast<double>(_repetitions);

      // Conservative decision rule: only accept "sorted wins" once a clear majority of repetitions
      // agree by a clear margin (see _sortedWinMarginFraction and _requiredSortedWinRatio).
      if (winRatio >= _requiredSortedWinRatio) {
        highCount = mid;
        AutoPasLog(TRACE, "SortingThresholdBenchmark search {} layout={} n={}: sorted won {}/{} reps → high={}",
                   newton3, layout, mid, outcome, _repetitions, highCount);
      } else {
        lowCount = mid + 1;
        AutoPasLog(TRACE, "SortingThresholdBenchmark search {} layout={} n={}: sorted won only {}/{} reps → low={}",
                   newton3, layout, mid, outcome, _repetitions, lowCount);
      }
    }
    AutoPasLog(DEBUG, "SortingThresholdBenchmark {} layout={} threshold={}", newton3, layout, lowCount);
    return lowCount;
  }
};

}  // namespace autopas
