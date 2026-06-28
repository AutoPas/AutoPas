/**
 * @file SortingThresholdBenchmark.h
 * @date 27.06.2026
 * @author hmeyran
 */

#pragma once

#include <algorithm>
#include <array>

namespace autopas {

using size_t = std::size_t;

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
   * Whether runBenchmark() has already been called.
   * Checked by LogicHandler to avoid re-running the benchmark on every iteration.
   */
  bool _hasRun{false};

  SortingThresholdBenchmark() { _thresholds.fill(25); }

  /**
   * Return the threshold for the given sortingDirection.
   * @param sortingDirection Normalized direction vector; the number of zero components selects the index.
   * @return Particle count below which the unsorted path is used.
   */
  size_t getThreshold(const std::array<double, 3> &sortingDirection) const {
    const auto zeroCount = std::count(std::begin(sortingDirection), std::end(sortingDirection), 0.0);
    return _thresholds[static_cast<size_t>(zeroCount)];
  }

  /**
   * Stub: sets _hasRun and fills all thresholds with the default value.
   * Replace the body of this function with the actual micro-benchmark once the wiring is validated.
   * @tparam Functor_T
   * @tparam Particle_T
   * @param functor
   * @param sortingCutoff
   */
  template <class Functor_T, class Particle_T>
  void runBenchmark(Functor_T & /*functor*/, double /*sortingCutoff*/) {
    _thresholds.fill(25);
    _hasRun = true;
  }

 private:
  /**
   * Per-direction-type threshold values.
   * Index = number of zero components in sortingDirection (0=Corner, 1=Edge, 2=Face).
   */
  std::array<size_t, 3> _thresholds{};
};

}  // namespace autopas
