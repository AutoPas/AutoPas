/**
 * @file SortingThresholdBenchmark.h
 * @date 27.06.2026
 * @author hmeyran
 */

#pragma once

#include <array>
#include <string_view>
#include <utility>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/logging/Logger.h"
#include "autopasTools/generators/UniformGenerator.h"

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
 */
class SortingThresholdBenchmark {
 public:
  /**
   * Initializes all thresholds to the compile-time default.
   */
  SortingThresholdBenchmark() { _thresholds.fill(25); }

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
   * Per-direction-type threshold values.
   * Index = number of zero components in sortingDirection (0=Corner, 1=Edge, 2=Face).
   */
  std::array<size_t, 3> _thresholds{};

  /**
   * Number of timed calls per repetition; amortizes timer overhead within a single rep.
   * @todo: Adjust this based on input size.
   */
  size_t _iterations = 100;

  /**
   * Number of independent measurement repetitions per particle count; the mean is taken over these.
   * @todo: Try out different numbers of repetitions.
   */
  size_t _repetitions = 25;
  /**
   * Upper bound on the particle count searched by the binary search.
   */
  size_t _maxParticles = 500;

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
   * @return {sorted_time_ns, unsorted_time_ns}, each averaged over the configured number of repetitions.
   */
  template <class Functor_T, class Particle_T>
  std::pair<size_t, size_t> executeRun(Functor_T &functor, size_t layout, size_t numParticles) {
    using BenchCell = FullParticleCell<Particle_T>;
    using BenchCF = internal::CellFunctor<BenchCell, Functor_T, false>;

    const Particle_T defaultParticle({0, 0, 0}, {0, 0, 0}, 0);
    const double cutoff = functor.getCutoff();
    const double invSqrt3 = 1. / sqrt(3.);
    const double invSqrt2 = 1. / sqrt(2.);
    BenchCF cellFunctor{functor, functor.getCutoff(), DataLayoutOption::soa, false};
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
    constexpr std::array<std::string_view, 3> layoutNames = {"Corner", "Edge", "Face"};
    utils::Timer sortedTimer, unsortedTimer;
    for (size_t i = 0; i < _repetitions; i++) {
      cell1.clear();
      cell2.clear();

      autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, cell1Low, cell1High,
                                                                    numParticles);
      autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, cell2Low, cell2High,
                                                                    numParticles);
      functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
      functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);

      long beforeUnsorted = unsortedTimer.getTotalTime();
      unsortedTimer.start();
      for (size_t j = 0; j < _iterations; j++) {
        // A sorting direction of (0, 0, 0) disables sorting.
        cellFunctor.processCellPair(cell1, cell2, {0., 0., 0.});
      }
      unsortedTimer.stop();

      long beforeSorted = sortedTimer.getTotalTime();
      sortedTimer.start();
      for (size_t j = 0; j < _iterations; j++) {
        cellFunctor.processCellPair(cell1, cell2, sortingDirection);
      }
      sortedTimer.stop();

      AutoPasLog(DEBUG, "SortingThresholdBenchmark rep {}/{} layout={} n={}: unsorted={}ns sorted={}ns", i + 1,
                 _repetitions, layoutNames[layout], numParticles, unsortedTimer.getTotalTime() - beforeUnsorted,
                 sortedTimer.getTotalTime() - beforeSorted);
    }

    const long meanSorted = sortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    const long meanUnsorted = unsortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    AutoPasLog(INFO, "SortingThresholdBenchmark layout={} n={}: mean unsorted={}ns mean sorted={}ns",
               layoutNames[layout], numParticles, meanUnsorted, meanSorted);
    return {static_cast<size_t>(meanSorted), static_cast<size_t>(meanUnsorted)};
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
    constexpr std::array<std::string_view, 3> layoutNames = {"Corner", "Edge", "Face"};
    size_t lowCount = 0;
    size_t highCount = _maxParticles;

    while (lowCount < highCount) {
      size_t mid = lowCount + (highCount - lowCount) / 2;

      auto [sortedT, unsortedT] = executeRun<Functor_T, Particle_T>(functor, layout, mid);
      if (sortedT < unsortedT) {
        highCount = mid;
        AutoPasLog(DEBUG, "SortingThresholdBenchmark search layout={} n={}: sorted({}ns) < unsorted({}ns) → high={}",
                   layoutNames[layout], mid, sortedT, unsortedT, highCount);
      } else {
        lowCount = mid + 1;
        AutoPasLog(DEBUG, "SortingThresholdBenchmark search layout={} n={}: sorted({}ns) >= unsorted({}ns) → low={}",
                   layoutNames[layout], mid, sortedT, unsortedT, lowCount);
      }
    }
    AutoPasLog(INFO, "SortingThresholdBenchmark layout={} threshold={}", layoutNames[layout], lowCount);
    return lowCount;
  }
};

}  // namespace autopas
