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
     * ExpertSel: uses runtime performance info and fuzzy logic with expert rules to select scheduling algorithms.
     */
    auto4omp_expertsel,

#ifdef AUTOPAS_USE_LB4OMP
    // LB4OMP's scheduling techniques:

    // Dynamic, non-adaptive OpenMP non-standard techniques:
    lb4omp_trapezoid_self_scheduling, // TSS

    // Dynamic, non-adaptive LB4OMP techniques:
    lb4omp_fixed_size_chunk, // FSC
    lb4omp_factoring, // FAC
    lb4omp_improved_factoring, // mFAC
    lb4omp_practical_factoring, // FAC2
    lb4omp_practical_weighted_factoring, // WF2
    lb4omp_tapering, // TAP
    lb4omp_modified_fixed_size_chunk, // mFSC
    lb4omp_trapezoid_factoring_self_scheduling, // TFSS
    lb4omp_fixedIncrease_self_scheduling, // FISS
    lb4omp_variable_increase_self_scheduling, // FISS
    lb4omp_random, // RND

    // Dynamic, adaptive LB4OMP techniques:
    lb4omp_bold, // BOLD
    lb4omp_adaptive_weighted_factoring, // AWF
    lb4omp_adaptive_weighted_factoring_B, // AWF-B
    lb4omp_adaptive_weighted_factoring_C, // AWF-C
    lb4omp_adaptive_weighted_factoring_D, // AWF-D
    lb4omp_adaptive_weighted_factoring_E, // AWF-E
    lb4omp_adaptive_factoring, // AF
    lb4omp_improved_adaptive_factoring, // mAF

    // Performance measurement for FSC, FAC, TAP, BOLD.
    lb4omp_profiling,
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
        {OpenMPKindOption::auto4omp_expertsel, "expertSel"},

#ifdef AUTOPAS_USE_LB4OMP
        // LB4OMP's scheduling techniques:
        {OpenMPKindOption::lb4omp_trapezoid_self_scheduling, "Trapezoid self scheduling (TSS)"},
        {OpenMPKindOption::lb4omp_fixed_size_chunk, "Fixed size chunk (FSC)"},
        {OpenMPKindOption::lb4omp_factoring, "Factoring (FAC)"},
        {OpenMPKindOption::lb4omp_improved_factoring, "Improved implementation of Factoring (mFAC)"},
        {OpenMPKindOption::lb4omp_practical_factoring, "Practical variant of factoring (FAC2)"},
        {OpenMPKindOption::lb4omp_practical_weighted_factoring, "Practical variant of weighted factoring (WF2)"},
        {OpenMPKindOption::lb4omp_tapering, "Tapering (TAP)"},
        {OpenMPKindOption::lb4omp_modified_fixed_size_chunk, "Modified Fixed size chunk (mFSC)"},
        {OpenMPKindOption::lb4omp_trapezoid_factoring_self_scheduling, "Trapezoid factoring self scheduling "
         "(TFSS)"},
        {OpenMPKindOption::lb4omp_fixedIncrease_self_scheduling, "Fixed increase self scheduling (FISS)"},
        {OpenMPKindOption::lb4omp_variable_increase_self_scheduling, "Variable increase self scheduling (FISS)"},
        {OpenMPKindOption::lb4omp_random, "Random (RND)"},
        {OpenMPKindOption::lb4omp_bold, "(BOLD)"},
        {OpenMPKindOption::lb4omp_adaptive_weighted_factoring, "Adaptive weighted factoring (AWF) "
         "for time-stepping applications"},
        {OpenMPKindOption::lb4omp_adaptive_weighted_factoring_B, "Variant B of adaptive weighted factoring "
         "(AWF-B)"},
        {OpenMPKindOption::lb4omp_adaptive_weighted_factoring_C, "Variant C of adaptive weighted factoring "
         "(AWF-C)"},
        {OpenMPKindOption::lb4omp_adaptive_weighted_factoring_D, "Variant D of adaptive weighted factoring "
         "(AWF-D)"},
        {OpenMPKindOption::lb4omp_adaptive_weighted_factoring_E, "Variant E of adaptive weighted factoring "
         "(AWF-E)"},
        {OpenMPKindOption::lb4omp_adaptive_factoring, "Adaptive factoring (AF)"},
        {OpenMPKindOption::lb4omp_improved_adaptive_factoring, "Improved implementation of Adaptive factoring "
         "(mAF)"},
        {OpenMPKindOption::lb4omp_profiling, ""},
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
   * OpenMP chunk size getter.
   * @return the current OpenMP chunk size
   */
  [[maybe_unused]] [[nodiscard]] int getChunkSize() const;

  /**
   * OpenMP chunk size getter for setting OpenMP's scheduling runtime variables.
   * @return the current OpenMP chunk size
   */
  [[maybe_unused]] [[nodiscard]] int getOMPChunkSize() const;

  /**
   * OpenMP chunk size setter.
   * @param s the new chunk size to use
   */
  [[maybe_unused]] void setChunkSize(int chunkSize);

  /**
   * OpenMP scheduling kind getter.
   * @return the current OpenMP scheduling kind, directly usable as an argument for OpenMP's schedule-clause
   */
  [[maybe_unused]] [[nodiscard]] OpenMPKindOption getKind() const;

  /**
   * OpenMP standard scheduling kind getter.
   * @return the current OpenMP chunk size
   */
  [[maybe_unused]] [[nodiscard]] omp_sched_t getOMPKind() const;

  /**
   * OpenMP scheduling kind setter.
   * @param k the new scheduling kind to use
   */
  [[maybe_unused]] void setKind(OpenMPKindOption kind);

  /**
   * Sets OpenMP's scheduling runtime variables.
   */
  [[maybe_unused]] inline void setSchedule() const { autopas_set_schedule(getOMPKind(), getOMPChunkSize()); }

  /**
   * Sets OpenMP's scheduling runtime variables.
   * @param k the new scheduling kind
   * @param s the new chunk size
   */
  [[maybe_unused]] inline void setSchedule(OpenMPKindOption kind, int chunkSize) {
    setKind(kind);
    setChunkSize(chunkSize);
    setSchedule();
  }

  /**
   * Tells whether the scheduling chunk size should be overwritten.
   * @return whether the scheduling chunk size should be overwritten
   */
  [[maybe_unused]] [[nodiscard]] bool overrideChunkSize() const;

  /**
   * Tells whether the scheduling kind is a manual LB4OMP scheduling technique.
   * @return whether the scheduling kind is a manual LB4OMP scheduling technique
   */
  [[maybe_unused]] [[nodiscard]] bool manualSchedulingTechnique() const;
};
}  // namespace autopas