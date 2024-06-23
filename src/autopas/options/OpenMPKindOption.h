/**
 * @file OpenMPKindOption.h
 * @author MehdiHachicha
 * @date 22.06.2024
 */

#pragma once

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {

/**
 * Tests if string s contains sub-string sub.
 * @param s the string
 * @param sub the sub-string
 * @return whether s contains sub
 */
static bool contains(const std::string &s, const std::string &sub) { return s.find(sub) != std::string::npos; }

/**
 * List of valid scheduling kind abbreviations.
 */
static const std::array<std::string, 36> validNameAbbreviations{
    // clang-format off
    "dyn",  "SS",   "guid", "GSS",  "run", "sta", "STATIC", "rand", "exh",   "bin",   "exp",
#ifdef AUTOPAS_USE_LB4OMP
    "prof", "fsc",  "FSC",  "tap",  "TAP", "fac", "FAC",    "bold", "BOLD",  "wf",    "WF",  "tfss", "TFSS",
    "fiss", "FISS", "viss", "VISS", "rnd", "RND", "trap",   "TSS",  "steal", "Steal", "af",  "AF"
#endif
    // clang-format on
};

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
    // LB4OMP's scheduling techniques [1, 2] (beware, technique names in the papers and Git README are outdated):
    /**
     * Profiling
     * Pass a path for storing the data with KMP_PROFILE_DATA.
     */
    lb4omp_profiling,

    /**
     * Fixed Size Chunk
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_fsc,

    /**
     * Modified Fixed Size Chunk
     */
    lb4omp_mfsc,

    /**
     * Tapering
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_tap,

    /**
     * Factoring
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_fac,

    /**
     * Improved Factoring
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_faca,

    /**
     * Bold
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_bold,

    /**
     * Practical Factoring
     */
    lb4omp_fac2,

    /**
     * Weighted Factoring
     */
    lb4omp_wf,

    /**
     * Adaptive Factoring
     */
    lb4omp_af,

    /**
     * Adaptive Weighted Factoring
     */
    lb4omp_awf,

    /**
     * Trapezoid Factoring Self Scheduling
     */
    lb4omp_tfss,

    /**
     * Fixed Increase Self Scheduling
     */
    lb4omp_fiss,

    /**
     * Variable Increase Self Scheduling
     */
    lb4omp_viss,

    /**
     * Random
     */
    lb4omp_rnd,

    // LB4OMP's scheduling techniques used by Auto4OMP (in addition to the standard scheduling kinds) [1, 2]:
    /**
     * Trapezoid Self Scheduling (from standard OpenMP)
     */
    lb4omp_trapezoidal,

    /**
     * Improved Practical Factoring
     */
    lb4omp_fac2a,

    /**
     * Static with Steal enabled (from standard OpenMP)
     */
    lb4omp_static_steal,

    /**
     * Adaptive Weighted Factoring Variant B
     */
    lb4omp_awf_b,

    /**
     * Adaptive Weighted Factoring Variant C
     */
    lb4omp_awf_c,

    /**
     * Adaptive Weighted Factoring Variant D
     */
    lb4omp_awf_d,

    /**
     * Adaptive Weighted Factoring Variant E
     */
    lb4omp_awf_e,

    /**
     * Improved Adaptive Factoring
     */
    lb4omp_af_a,
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
        {OpenMPKindOption::lb4omp_trapezoidal, "trapezoidal"},    // Trapezoid Self Scheduling
        {OpenMPKindOption::lb4omp_fac2a, "fac2a"},                // Improved Practical Factoring
        {OpenMPKindOption::lb4omp_static_steal, "static_steal"},  // Static with Steal enabled
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

  /**
   * Tells whether a given kind name is valid.
   * @param name the name to test
   * @return whether the name is a valid scheduling kind option
   */
  static bool valid(const std::string &name) {
    return std::any_of(validNameAbbreviations.begin(), validNameAbbreviations.end(),
                       [=](const std::string &abbreviation) { return contains(name, abbreviation); });
  }

  /**
   * Parse the scheduling kind option from a string.
   * @param the name of the scheduling kind
   * @return the corresponding AutoPas option
   */
  static OpenMPKindOption parse(const std::string &name) {
    if (contains(name, "dyn") || contains(name, "SS"))
      return OpenMPKindOption::omp_dynamic;
    else if (contains(name, "guid") || contains(name, "GSS"))
      return OpenMPKindOption::omp_guided;
    else if (contains(name, "run"))
      return OpenMPKindOption::omp_runtime;
    else if (contains(name, "sta") || contains(name, "STATIC"))
      return OpenMPKindOption::omp_static;
    else if (contains(name, "rand"))
      return OpenMPKindOption::auto4omp_randomsel;
    else if (contains(name, "exh"))
      return OpenMPKindOption::auto4omp_exhaustivesel;
    else if (contains(name, "bin"))
      return OpenMPKindOption::auto4omp_binarySearch;
    else if (contains(name, "exp"))
      return OpenMPKindOption::auto4omp_expertsel;
#ifdef AUTOPAS_USE_LB4OMP  // LB4OMP's scheduling techniques.
    else if (contains(name, "prof"))
      return OpenMPKindOption::lb4omp_profiling;
    else if (contains(name, "mfsc") || contains(name, "mFSC"))
      return OpenMPKindOption::lb4omp_mfsc;
    else if (contains(name, "fsc") || contains(name, "FSC"))
      return OpenMPKindOption::lb4omp_fsc;
    else if (contains(name, "tap") || contains(name, "TAP"))
      return OpenMPKindOption::lb4omp_tap;
    else if (contains(name, "fac2a") || contains(name, "mFAC2"))
      return OpenMPKindOption::lb4omp_fac2a;
    else if (contains(name, "fac2") || contains(name, "FAC2"))
      return OpenMPKindOption::lb4omp_fac2;
    else if (contains(name, "faca") || contains(name, "mFAC"))
      return OpenMPKindOption::lb4omp_faca;
    else if (contains(name, "fac") || contains(name, "FAC"))
      return OpenMPKindOption::lb4omp_fac;
    else if (contains(name, "bold") || contains(name, "BOLD"))
      return OpenMPKindOption::lb4omp_bold;
    else if (contains(name, "awf_b") || contains(name, "AWF-B"))
      return OpenMPKindOption::lb4omp_awf_b;
    else if (contains(name, "awf_c") || contains(name, "AWF-C"))
      return OpenMPKindOption::lb4omp_awf_c;
    else if (contains(name, "awf_d") || contains(name, "AWF-D"))
      return OpenMPKindOption::lb4omp_awf_d;
    else if (contains(name, "awf_e") || contains(name, "AWF-E"))
      return OpenMPKindOption::lb4omp_awf_e;
    else if (contains(name, "awf") || contains(name, "AWF"))
      return OpenMPKindOption::lb4omp_awf;
    else if (contains(name, "wf") || contains(name, "WF"))
      return OpenMPKindOption::lb4omp_wf;
    else if (contains(name, "tfss") || contains(name, "TFSS"))
      return OpenMPKindOption::lb4omp_tfss;
    else if (contains(name, "fiss") || contains(name, "FISS"))
      return OpenMPKindOption::lb4omp_fiss;
    else if (contains(name, "viss") || contains(name, "VISS"))
      return OpenMPKindOption::lb4omp_viss;
    else if (contains(name, "rnd") || contains(name, "RND"))
      return OpenMPKindOption::lb4omp_rnd;
    else if (contains(name, "trap") || contains(name, "TSS"))
      return OpenMPKindOption::lb4omp_trapezoidal;
    else if (contains(name, "steal") || contains(name, "Steal"))
      return OpenMPKindOption::lb4omp_static_steal;
    else if (contains(name, "af_a") || contains(name, "mAF"))
      return OpenMPKindOption::lb4omp_af_a;
    else if (contains(name, "af") || contains(name, "AF"))
      return OpenMPKindOption::lb4omp_af;
#endif
    else  // Default.
      return OpenMPKindOption::omp_auto;
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas

/*
 * Sources:
 * [1] https://www.computer.org/csdl/journal/td/2022/04/09524500/1wpqIcNI6YM
 * [2] https://ieeexplore.ieee.org/document/9825675
 */
