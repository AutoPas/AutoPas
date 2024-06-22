/**
 * @file OpenMPConfigurator.h
 * @author MehdiHachicha
 * @date 12.03.2024
 */

#pragma once

#include <cstddef>
#include <set>

#include "autopas/options/Option.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * Class representing OpenMP's scheduling kind choices.
 */
class OpenMPKindOption : public Option<OpenMPKindOption> {
 public:
  /**
   * Possible choices for OpenMP's scheduling kind.
   */
  enum Value {
    // Standard OpenMP's scheduling kinds:

    /**
     * auto: automatic scheduling decisions.
     */
    omp_auto,

    /**
     * dynamic: iterations are distributed to the threads in chunks.
     * When a thread finishes its chunk, it requests a new one.
     * The chunk size remains constant.
     */
    omp_dynamic,

    /**
     * guided: iterations are distributed to the threads in chunks.
     * When a thread finishes its chunk, it requests a new one.
     * The chunk size starts large, and decreases over time towards the minimum set by the chunk size argument.
     */
    omp_guided,

    /**
     * runtime: uses the scheduling kind set by the OMP_SCHEDULE environment variable.
     */
    omp_runtime,

    /**
     * omp_static: iterations are distributed to the threads in chunks.
     * The chunks are all distributed in advance (round-robin).
     * If iterations runtimes differ significantly, some threads may finish too early.
     * In such cases, use dynamic scheduling instead.
     */
    omp_static,

    // Auto4OMP's automated selection methods:

    /**
     * RandomSel: may reselect a random scheduling algorithm at the end of each time-step.
     * The probability of reselection increases with greater load imbalance.
     */
    auto4omp_randomsel,

    /**
     * ExhaustiveSel: Tries each algorithm once, then selects the fastest one and keeps it for a while.
     * Load imbalance re-triggers selection.
     */
    auto4omp_exhaustivesel,

    /**
     * Binary Search
     */
    auto4omp_binarySearch,

    /**
     * ExpertSel: uses runtime performance info and fuzzy logic with expert rules to select scheduling algorithms.
     */
    auto4omp_expertsel,

#ifdef AUTOPAS_USE_LB4OMP
    // LB4OMP's scheduling techniques (beware, technique names in LB4OMP's README are outdated) [1, 3]:
    lb4omp_profiling,  // Profiling                    // Point KMP_PROFILE_DATA to where to store the data
    lb4omp_fsc,        // Fixed Size Chunk             // Requires profiling, set KMP_PROFILE_DATA
    lb4omp_mfsc,       // Modified Fixed Size Chunk
    lb4omp_tap,        // Tapering                     // Requires profiling, set KMP_PROFILE_DATA
    lb4omp_fac,        // Factoring                    // Requires profiling, set KMP_PROFILE_DATA
    lb4omp_faca,       // Improved Factoring           // Requires profiling, set KMP_PROFILE_DATA
    lb4omp_bold,       // Bold                         // Requires profiling, set KMP_PROFILE_DATA
    lb4omp_fac2,       // Practical Factoring
    lb4omp_wf,         // Weighted Factoring
    lb4omp_af,         // Adaptive Factoring
    lb4omp_awf,        // Adaptive Weighted Factoring
    lb4omp_tfss,       // Trapezoid Factoring Self Scheduling
    lb4omp_fiss,       // Fixed Increase Self Scheduling
    lb4omp_viss,       // Variable Increase Self Scheduling
    lb4omp_rnd,        // Random

    // LB4OMP's scheduling techniques used by Auto4OMP (in addition to the standard scheduling kinds) [1, 2, 3]:
    lb4omp_trapezoidal,   // Trapezoid Self Scheduling (from standard OpenMP)
    lb4omp_fac2a,         // Improved Practical Factoring
    lb4omp_static_steal,  // Static with Steal enabled (from standard OpenMP)
    lb4omp_awf_b,         // Adaptive Weighted Factoring Variant B
    lb4omp_awf_c,         // Adaptive Weighted Factoring Variant C
    lb4omp_awf_d,         // Adaptive Weighted Factoring Variant D
    lb4omp_awf_e,         // Adaptive Weighted Factoring Variant E
    lb4omp_af_a,          // Improved Adaptive Factoring
#endif
  };

  /**
   * Constructor.
   */
  OpenMPKindOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  // NOLINTNEXTLINE: marking explicit triggers errors elsewhere. Keep non-explicit for now.
  constexpr OpenMPKindOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  // NOLINTNEXTLINE: marking explicit triggers errors elsewhere. Keep non-explicit for now.
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<OpenMPKindOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of OpenMPKindOption.
   * @return map option -> string representation
   */
  static std::map<OpenMPKindOption, std::string> getOptionNames() {
    return {
        // Standard OpenMP's scheduling kinds:
        {OpenMPKindOption::omp_auto, "auto"},
        {OpenMPKindOption::omp_dynamic, "dynamic"},
        {OpenMPKindOption::omp_guided, "guided"},
        {OpenMPKindOption::omp_runtime, "runtime"},
        {OpenMPKindOption::omp_static, "static"},

        // Auto4OMP's automated selection methods:
        {OpenMPKindOption::auto4omp_randomsel, "randomSel"},
        {OpenMPKindOption::auto4omp_exhaustivesel, "exhaustiveSel"},
        {OpenMPKindOption::auto4omp_binarySearch, "binarySearch"},
        {OpenMPKindOption::auto4omp_expertsel, "expertSel"},

#ifdef AUTOPAS_USE_LB4OMP
        // LB4OMP's scheduling techniques (beware, technique names in LB4OMP's README are outdated):
        {OpenMPKindOption::lb4omp_profiling, "profiling"},  // Profiling
        {OpenMPKindOption::lb4omp_fsc, "fsc"},              // Fixed Size Chunk             // Requires profiling
        {OpenMPKindOption::lb4omp_mfsc, "mfsc"},            // Modified Fixed Size Chunk
        {OpenMPKindOption::lb4omp_tap, "tap"},              // Tapering                     // Requires profiling
        {OpenMPKindOption::lb4omp_fac, "fac"},              // Factoring                    // Requires profiling
        {OpenMPKindOption::lb4omp_faca, "faca"},            // Improved Factoring           // Requires profiling
        {OpenMPKindOption::lb4omp_bold, "bold"},            // Bold                         // Requires profiling
        {OpenMPKindOption::lb4omp_fac2, "fac2"},            // Practical Factoring
        {OpenMPKindOption::lb4omp_wf, "wf"},                // Weighted Factoring
        {OpenMPKindOption::lb4omp_af, "af"},                // Adaptive Factoring
        {OpenMPKindOption::lb4omp_awf, "awf"},              // Adaptive Weighted Factoring
        {OpenMPKindOption::lb4omp_tfss, "tfss"},            // Trapezoid Factoring Self Scheduling
        {OpenMPKindOption::lb4omp_fiss, "fiss"},            // Fixed Increase Self Scheduling
        {OpenMPKindOption::lb4omp_viss, "viss"},            // Variable Increase Self Scheduling
        {OpenMPKindOption::lb4omp_rnd, "rnd"},              // Random

        // LB4OMP's scheduling techniques used by Auto4OMP (in addition to the standard scheduling kinds):
        {OpenMPKindOption::lb4omp_trapezoidal, "trapezoidal"},    // Trapezoid Self Scheduling (from standard OpenMP)
        {OpenMPKindOption::lb4omp_fac2a, "fac2a"},                // Improved Practical Factoring
        {OpenMPKindOption::lb4omp_static_steal, "static_steal"},  // Static with Steal enabled (from standard OpenMP)
        {OpenMPKindOption::lb4omp_awf_b, "awf_b"},                // Adaptive Weighted Factoring Variant B
        {OpenMPKindOption::lb4omp_awf_c, "awf_c"},                // Adaptive Weighted Factoring Variant C
        {OpenMPKindOption::lb4omp_awf_d, "awf_d"},                // Adaptive Weighted Factoring Variant D
        {OpenMPKindOption::lb4omp_awf_e, "awf_e"},                // Adaptive Weighted Factoring Variant E
        {OpenMPKindOption::lb4omp_af_a, "af_a"},                  // Improved Adaptive Factoring
#endif
    };
  };

  /**
   * Checks if a given kind is in a given list of kinds.
   * @return whether the given kind is in the given list of kinds
   */
  template <typename... KindList>
  static bool in(OpenMPKindOption k, KindList... list) {
    return (... || (k == list));
  }

 private:
  Value _value{Value(-1)};
};

/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-chunk-size <int>
 */
extern int openMPDefaultChunkSize;

/**
 * OpenMP default scheduling kind for manual testing.
 * md-flexible: set via command-line option --openmp-kind <OpenMPKindOption>
 */
extern OpenMPKindOption openMPDefaultKind;

/**
 * This class provides configurable parameters for OpenMP.
 */
class OpenMPConfigurator {
 private:
  OpenMPKindOption _kind = openMPDefaultKind;
  int _chunkSize = openMPDefaultChunkSize;

 public:
  /**
   * OpenMP configurator default constructor.
   */
  [[maybe_unused]] OpenMPConfigurator();

  /**
   * OpenMP configurator constructor.
   * @param s chunk size used in OpenMP's loop scheduling
   */
  [[maybe_unused]] explicit OpenMPConfigurator(OpenMPKindOption kind, int chunkSize);

  /**
   * AutoPas OpenMP configurator chunk size getter.
   * @return the current OpenMP chunk size
   */
  [[maybe_unused]] [[nodiscard]] int getChunkSize() const;

  /**
   * OpenMP chunk size getter for setting OpenMP's scheduling runtime variables.
   * @return the current OpenMP chunk size, directly usable in OpenMP's schedule setter
   */
  [[maybe_unused]] [[nodiscard]] int getOMPChunkSize() const;

  /**
   * AutoPas OpenMP configurator chunk size setter.
   * @param chunkSize the new chunk size to use
   */
  [[maybe_unused]] void setChunkSize(int chunkSize);

  /**
   * AutoPas OpenMP configurator scheduling kind getter.
   * @return the current OpenMP scheduling kind
   */
  [[maybe_unused]] [[nodiscard]] OpenMPKindOption getKind() const;

  /**
   * OpenMP scheduling kind getter for setting OpenMP's scheduling runtime variables.
   * @return the current OpenMP kind, directly usable in OpenMP's schedule setter
   */
  [[maybe_unused]] [[nodiscard]] omp_sched_t getOMPKind() const;

  /**
   * AutoPas OpenMP configurator scheduling kind setter.
   * @param kind the new scheduling kind to use
   */
  [[maybe_unused]] void setKind(OpenMPKindOption kind);

  /**
   * Tells whether the scheduling chunk size should be overwritten.
   * @return whether the scheduling chunk size should be overwritten
   */
  [[maybe_unused]] [[nodiscard]] bool overrideChunkSize() const;

  /**
   * Tells whether the scheduling kind is a standard OpenMP kind.
   * @return whether the scheduling kind is a standard OpenMP kind
   */
  [[maybe_unused]] [[nodiscard]] bool standard() const;
};  // class OpenMPConfigurator

/**
 * Sets OpenMP's runtime schedule from a given OpenMP configurator.
 * schedule(runtime) will then use them for the traversal in the concerned calling thread.
 * @param ompConfig the OpenMP configurator
 */
inline void autopas_set_schedule(autopas::OpenMPConfigurator ompConfig) {
  // If the configurator is set to omp_runtime, users are assumed to have set OMP_SCHEDULE manually.
  if (ompConfig.getKind() == OpenMPKindOption::omp_runtime) return;

  if (ompConfig.standard()) {
    autopas_set_schedule(ompConfig.getOMPKind(), ompConfig.getOMPChunkSize());
  } else {
    autopas_auto4omp_set_schedule(ompConfig.getOMPKind(), ompConfig.getOMPChunkSize());
  }
}  // void autopas_set_schedule
}  // namespace autopas

/*
 * Sources:
 * [1] https://www.computer.org/csdl/journal/td/2022/04/09524500/1wpqIcNI6YM
 * [2] https://ieeexplore.ieee.org/document/9825675
 */
