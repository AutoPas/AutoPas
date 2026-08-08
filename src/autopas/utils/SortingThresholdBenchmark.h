/**
 * @file SortingThresholdBenchmark.h
 * @date 27.06.2026
 * @author hmeyran
 */

#pragma once

#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/options/SortingDirectionOption.h"
#include "autopas/utils/SortingThresholdInfo2B.h"
#include "autopas/utils/SortingThresholdInfoSingle.h"
#include "autopas/utils/Timer.h"
#include "autopas/utils/generators/UniformGenerator.h"
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
 * The search performed in runSearch() is deliberately biased towards a conservative (higher) threshold, to avoid noisy
 * measurements influencing the threshold to be too low. A too low threshold would actively regress performance while a
 * too high threshold would at worst still perform the same as baseline. One goal of these thresholds is to avoid any
 * performance regression against the non sorted case.
 *
 * DEBUG BRANCH (debug/auto-threshold-side-effect). This class carries two environment-variable hooks that exist
 * only to split the md-flexible spinodal-decomposition regression at low cellSizeFactor into its two possible
 * causes. In that scenario ForceUpdateTotal is 10% higher under the calibrated thresholds than under any forced
 * threshold from 77 to 10000, even though the calibrated Newton3-on values (corner 77, edge 79, face 96) make the
 * calibrated run a behavioural subset of the forced-77 run, which shows no regression at all. Either the six values
 * are somehow responsible after all, or merely running the calibration perturbs the simulation that follows it.
 * The hooks isolate one from the other:
 *
 *   AUTOPAS_DEBUG_SORTING_THRESHOLD_OVERRIDE_SOA=<n>
 *   AUTOPAS_DEBUG_SORTING_THRESHOLD_OVERRIDE_AOS=<n>
 *     Run the calibration as usual, then discard its result and store a uniform <n> instead. Keeps every side
 *     effect of having calibrated, removes every decision the calibration would have influenced. A regression that
 *     survives this is not caused by the threshold values.
 *
 *   AUTOPAS_DEBUG_SORTING_THRESHOLD_INJECT_SOA=noN3Face,noN3Edge,noN3Corner,n3Face,n3Edge,n3Corner
 *   AUTOPAS_DEBUG_SORTING_THRESHOLD_INJECT_AOS=noN3Face,noN3Edge,noN3Corner,n3Face,n3Edge,n3Corner
 *     Install the six per-(Newton3, direction) values directly and mark the corresponding benchmark half as
 *     already run, so it is skipped. Keeps the values, removes the side effect. A regression that appears here is
 *     caused by the threshold values.
 *
 * Both halves are independent, and an unset or malformed variable leaves stock behaviour untouched. The applied
 * thresholds are logged at INFO in every case, so a run's own output records which values it used.
 */
class SortingThresholdBenchmark {
 public:
  /**
   * Constructs a benchmark whose thresholds default to the given uniform values until/unless a benchmark run
   * overwrites them. The defaults are stored as a SortingThresholdInfoSingle rather than a SortingThresholdInfo2B so
   * that getAoSThreshold()/getSoAThreshold() are always safe to hand to either CellFunctor or CellFunctor3B.
   * @param aosSortingThresholdDefault Uniform AoS threshold used before/without a benchmark run.
   * @param soaSortingThresholdDefault Uniform SoA threshold used before/without a benchmark run.
   */
  SortingThresholdBenchmark(size_t aosSortingThresholdDefault, size_t soaSortingThresholdDefault)
      : _soaThresholds(std::make_shared<const SortingThresholdInfoSingle>(soaSortingThresholdDefault)),
        _aosThresholds(std::make_shared<const SortingThresholdInfoSingle>(aosSortingThresholdDefault)) {
    applyDebugInjection();
  }

  /**
   * Return all per-Newton3-state, per-direction-type SoA sorting thresholds if hasRunSoA() is true, otherwise the
   * uniform default passed to the constructor.
   * @return Shared pointer to the internal threshold values.
   */
  [[nodiscard]] std::shared_ptr<const SortingThresholdInfoInterface> getSoAThreshold() const { return _soaThresholds; }

  /**
   * Return all per-Newton3-state, per-direction-type AoS pair-sorting thresholds if hasRunAoS() is true, otherwise
   * the uniform default passed to the constructor.
   * @return Shared pointer to the internal threshold values.
   */
  [[nodiscard]] std::shared_ptr<const SortingThresholdInfoInterface> getAoSThreshold() const { return _aosThresholds; }

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
   * @param defaultParticle Template particle (e.g. sampled from the particle container) whose non-positional
   * properties are copied onto every particle generated for the benchmark cells.
   */
  template <class Functor_T, class Particle_T>
  void runSoABenchmark(Functor_T &functor, const Particle_T &defaultParticle) {
    SortingThresholdInfo2B thresholds = SortingThresholdInfo2B(0);
    for (const auto &n3 : Newton3Option::getAllOptions()) {
      for (const auto &cellDirection : SortingDirectionOption::getAllOptions()) {
        thresholds.setThresholdByOption(
            n3, cellDirection, runSearch<Functor_T, Particle_T, true>(functor, defaultParticle, cellDirection, n3));
      }
    }
    _soaThresholds = std::make_shared<const SortingThresholdInfo2B>(thresholds);
    _hasRunSoA = true;
    logThresholds("SoA", "calibrated", thresholds);
    applyDebugOverride(_kEnvOverrideSoA, "SoA", _soaThresholds);
  }

  /**
   * Runs the micro-benchmark for the AoS sorting threshold, sweeping both Newton3 states and all three direction
   * types, and storing the resulting thresholds. Sets hasRunAoS() to true.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @param functor Functor instance used to drive the benchmark cells.
   * @param defaultParticle Template particle (e.g. sampled from the particle container) whose non-positional
   * properties are copied onto every particle generated for the benchmark cells.
   */
  template <class Functor_T, class Particle_T>
  void runAoSBenchmark(Functor_T &functor, const Particle_T &defaultParticle) {
    SortingThresholdInfo2B thresholds = SortingThresholdInfo2B(0);
    for (const auto &n3 : Newton3Option::getAllOptions()) {
      for (const auto &cellDirection : SortingDirectionOption::getAllOptions()) {
        thresholds.setThresholdByOption(
            n3, cellDirection, runSearch<Functor_T, Particle_T, false>(functor, defaultParticle, cellDirection, n3));
      }
    }
    _aosThresholds = std::make_shared<const SortingThresholdInfo2B>(thresholds);
    _hasRunAoS = true;
    logThresholds("AoS", "calibrated", thresholds);
    applyDebugOverride(_kEnvOverrideAoS, "AoS", _aosThresholds);
  }

 private:
  /**
   * Environment variable read by applyDebugOverride() for the SoA half. See the class documentation.
   */
  static constexpr const char *_kEnvOverrideSoA = "AUTOPAS_DEBUG_SORTING_THRESHOLD_OVERRIDE_SOA";
  /**
   * Environment variable read by applyDebugOverride() for the AoS half. See the class documentation.
   */
  static constexpr const char *_kEnvOverrideAoS = "AUTOPAS_DEBUG_SORTING_THRESHOLD_OVERRIDE_AOS";
  /**
   * Environment variable read by applyDebugInjection() for the SoA half. See the class documentation.
   */
  static constexpr const char *_kEnvInjectSoA = "AUTOPAS_DEBUG_SORTING_THRESHOLD_INJECT_SOA";
  /**
   * Environment variable read by applyDebugInjection() for the AoS half. See the class documentation.
   */
  static constexpr const char *_kEnvInjectAoS = "AUTOPAS_DEBUG_SORTING_THRESHOLD_INJECT_AOS";

  /**
   * Logs one set of per-(Newton3, direction) thresholds at INFO so a run records the values it used.
   * @param layout "SoA" or "AoS".
   * @param origin How the values were obtained, e.g. "calibrated" or "injected".
   * @param thresholds Values to log.
   */
  static void logThresholds(const char *layout, const char *origin, const SortingThresholdInfo2B &thresholds) {
    AutoPasLog(INFO, "SortingThresholdBenchmark {} thresholds ({}): noN3(face,edge,corner)=({},{},{}) n3=({},{},{})",
               layout, origin, thresholds.noN3FaceThreshold, thresholds.noN3EdgeThreshold,
               thresholds.noN3CornerThreshold, thresholds.n3FaceThreshold, thresholds.n3EdgeThreshold,
               thresholds.n3CornerThreshold);
  }

  /**
   * Debug hook: parse a single non-negative integer out of the environment.
   * @param name Environment variable to read.
   * @return The parsed value, or nullopt if unset or malformed.
   */
  static std::optional<size_t> parseDebugUniform(const char *name) {
    const char *raw = std::getenv(name);
    if (raw == nullptr) {
      return std::nullopt;
    }
    try {
      return static_cast<size_t>(std::stoull(raw));
    } catch (const std::exception &) {
      AutoPasLog(WARN, "SortingThresholdBenchmark: ignoring malformed {}='{}', expected one integer", name, raw);
      return std::nullopt;
    }
  }

  /**
   * Debug hook: parse six comma-separated non-negative integers out of the environment, in
   * SortingThresholdInfo2B constructor order (noN3Face, noN3Edge, noN3Corner, n3Face, n3Edge, n3Corner).
   * @param name Environment variable to read.
   * @return The parsed thresholds, or nullopt if unset or malformed.
   */
  static std::optional<SortingThresholdInfo2B> parseDebugPerConfig(const char *name) {
    const char *raw = std::getenv(name);
    if (raw == nullptr) {
      return std::nullopt;
    }
    std::vector<size_t> values;
    std::stringstream stream(raw);
    std::string token;
    while (std::getline(stream, token, ',')) {
      try {
        values.push_back(static_cast<size_t>(std::stoull(token)));
      } catch (const std::exception &) {
        values.clear();
        break;
      }
    }
    if (values.size() != 6) {
      AutoPasLog(WARN,
                 "SortingThresholdBenchmark: ignoring malformed {}='{}', expected six comma-separated integers "
                 "(noN3Face,noN3Edge,noN3Corner,n3Face,n3Edge,n3Corner)",
                 name, raw);
      return std::nullopt;
    }
    return SortingThresholdInfo2B(values[0], values[1], values[2], values[3], values[4], values[5]);
  }

  /**
   * Debug hook: if the injection variables are set, install their values and mark the corresponding benchmark half
   * as already run so runSoABenchmark()/runAoSBenchmark() are skipped. Called from the constructor, i.e. before any
   * calibration can happen. Isolates the effect of the calibrated values from the effect of calibrating.
   */
  void applyDebugInjection() {
    if (const auto injected = parseDebugPerConfig(_kEnvInjectSoA)) {
      _soaThresholds = std::make_shared<const SortingThresholdInfo2B>(*injected);
      _hasRunSoA = true;
      logThresholds("SoA", "injected, calibration skipped", *injected);
    }
    if (const auto injected = parseDebugPerConfig(_kEnvInjectAoS)) {
      _aosThresholds = std::make_shared<const SortingThresholdInfo2B>(*injected);
      _hasRunAoS = true;
      logThresholds("AoS", "injected, calibration skipped", *injected);
    }
  }

  /**
   * Debug hook: if the override variable is set, replace freshly calibrated thresholds with a uniform value.
   * Called at the end of a benchmark run, so the calibration itself still happened. Isolates the effect of
   * calibrating from the effect of the calibrated values.
   * @param envName Environment variable holding the uniform replacement value.
   * @param layout "SoA" or "AoS", for the log message.
   * @param thresholds Threshold member to overwrite.
   */
  static void applyDebugOverride(const char *envName, const char *layout,
                                 std::shared_ptr<const SortingThresholdInfoInterface> &thresholds) {
    if (const auto uniform = parseDebugUniform(envName)) {
      thresholds = std::make_shared<const SortingThresholdInfoSingle>(*uniform);
      AutoPasLog(INFO, "SortingThresholdBenchmark {} thresholds (calibrated result discarded via {}): uniform {}",
                 layout, envName, *uniform);
    }
  }

  /**
   * Set to true by runSoABenchmark() once the SoA benchmark has completed.
   */
  bool _hasRunSoA{false};

  /**
   * Set to true by runAoSBenchmark() once the AoS benchmark has completed.
   */
  bool _hasRunAoS{false};

  /**
   * SoA sorting threshold value. A SortingThresholdInfoSingle wrapping the constructor-provided default until
   * runSoABenchmark() replaces it with a SortingThresholdInfo2B.
   * Held as a shared_ptr so that getSoAThreshold() can hand out the same instance on every call once
   * runSoABenchmark() has completed, instead of forcing callers to reallocate on every access.
   */
  std::shared_ptr<const SortingThresholdInfoInterface> _soaThresholds;

  /**
   * AoS pair-sorting threshold value. A SortingThresholdInfoSingle wrapping the constructor-provided default until
   * runAoSBenchmark() replaces it with a SortingThresholdInfo2B.
   * Held as a shared_ptr so that getAoSThreshold() can hand out the same instance on every call once
   * runAoSBenchmark() has completed, instead of forcing callers to reallocate on every access.
   */
  std::shared_ptr<const SortingThresholdInfoInterface> _aosThresholds;

  /**
   * Number of timed calls per repetition to get stable measurement.
   */
  const size_t _iterations = 25;

  /**
   * Number of independent measurement repetitions per particle count.
   */
  const size_t _repetitions = 100;
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
   * Used to avoid noisy measurements influencing the data.
   */
  const double _sortedWinMarginFraction = 0.05;

  /**
   * Minimum fraction of repetitions within a single executeRun() call that must count as a clear win for
   * sorting (see _sortedWinMarginFraction) before runSearch() accepts "sorted is faster" for the tested
   * particle count.
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
   * The sorted path is controlled by passing a non zero sorting Direction.
   * @tparam Functor_T Pairwise functor type.
   * @tparam Particle_T Particle type.
   * @tparam useSoA Selects the benchmarked data layout: true measures the SoA path, false the AoS path.
   * @param functor Functor instance used to drive the benchmark.
   * @param defaultParticle Template particle whose non-positional properties are copied onto every particle
   * generated for the benchmark cells.
   * @param cellDirection Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param numParticles Combined number of particles placed in the two cells, with 40% of the particles landing in each
   * cell and 20% being randomly distributed between the cells.
   * @param newton3 Whether the benchmarked CellFunctor should apply Newton3.
   * @return The number of repetitions (out of _repetitions) in which the sorted path was measured as faster than
   * the unsorted path by at least _sortedWinMarginFraction.
   */
  template <class Functor_T, class Particle_T, bool useSoA>
  size_t executeRun(Functor_T &functor, const Particle_T &defaultParticle, SortingDirectionOption cellDirection,
                    size_t numParticles, Newton3Option newton3) {
    using BenchCell = FullParticleCell<Particle_T>;
    using BenchCF = internal::CellFunctor<BenchCell, Functor_T, true>;

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

    switch (cellDirection) {
      case SortingDirectionOption::corner:
        cell2Low = {cutoff, cutoff, cutoff};
        cell2High = {2. * cutoff, 2. * cutoff, 2. * cutoff};
        sortingDirection = {invSqrt3, invSqrt3, invSqrt3};
        break;
      case SortingDirectionOption::edge:
        cell2Low = {cutoff, cutoff, 0.};
        cell2High = {2. * cutoff, 2. * cutoff, cutoff};
        sortingDirection = {invSqrt2, invSqrt2, 0.};
        break;
      case SortingDirectionOption::face:
        cell2Low = {cutoff, 0., 0.};
        cell2High = {2 * cutoff, cutoff, cutoff};
        sortingDirection = {1, 0., 0.};
        break;
      default:
        utils::ExceptionHandler::exception("Cell Direction {} is not a valid/supported Direction!", cellDirection);
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

      AutoPasLog(TRACE, "SortingThresholdBenchmark rep {}/{} cell direction={} n={}: unsorted={}ns sorted={}ns", i + 1,
                 _repetitions, cellDirection, numParticles, unsortedDelta, sortedDelta);

      // A repetition only counts as a "sorted win" if it clears the margin: see _sortedWinMarginFraction.
      if (static_cast<double>(sortedDelta) < static_cast<double>(unsortedDelta) * (1. - _sortedWinMarginFraction)) {
        ++sortedWins;
      }
    }

    const long meanSorted = sortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    const long meanUnsorted = unsortedTimer.getTotalTime() / static_cast<long>(_repetitions);
    AutoPasLog(TRACE,
               "SortingThresholdBenchmark cell direction={} n={}: mean unsorted={}ns mean sorted={}ns sortedWins={}/{}",
               cellDirection, numParticles, meanUnsorted, meanSorted, sortedWins, _repetitions);
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
   * @param cellDirection Direction-type index (0=Corner, 1=Edge, 2=Face).
   * @param newton3 Whether the benchmarked CellFunctor should apply Newton3.
   * @return Smallest particle count at which sorted beats unsorted, or the upper search bound if never.
   */
  template <class Functor_T, class Particle_T, bool useSoA>
  size_t runSearch(Functor_T &functor, const Particle_T &defaultParticle, SortingDirectionOption cellDirection,
                   Newton3Option newton3) {
    size_t lowCount = 0;
    size_t highCount = useSoA ? _maxSoAParticles : _maxAoSParticles;

    while (lowCount < highCount) {
      size_t mid = lowCount + (highCount - lowCount) / 2;

      const auto outcome =
          executeRun<Functor_T, Particle_T, useSoA>(functor, defaultParticle, cellDirection, mid, newton3);
      const double winRatio = static_cast<double>(outcome) / static_cast<double>(_repetitions);

      // Conservative decision rule: only accept "sorted wins" once a clear majority of repetitions
      // agree by a clear margin (see _sortedWinMarginFraction and _requiredSortedWinRatio).
      if (winRatio >= _requiredSortedWinRatio) {
        highCount = mid;
        AutoPasLog(TRACE, "SortingThresholdBenchmark search {} cell direction={} n={}: sorted won {}/{} reps → high={}",
                   newton3, cellDirection, mid, outcome, _repetitions, highCount);
      } else {
        lowCount = mid + 1;
        AutoPasLog(TRACE,
                   "SortingThresholdBenchmark search {} cell direction={} n={}: sorted won only {}/{} reps → low={}",
                   newton3, cellDirection, mid, outcome, _repetitions, lowCount);
      }
    }
    AutoPasLog(DEBUG, "SortingThresholdBenchmark {} cell direction={} threshold={}", newton3, cellDirection, lowCount);
    return lowCount;
  }
};

}  // namespace autopas
