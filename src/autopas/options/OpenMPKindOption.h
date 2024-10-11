/**
 * @file OpenMPKindOption.h
 * @author MehdiHachicha
 * @date 22.06.2024
 */

#pragma once

#include "autopas/options/Option.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {
inline namespace options {

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
     * Standard OpenMP auto:
     * defaults to an analytical variant of guided.
     */
    omp_auto,

    /**
     * Standard OpenMP dynamic: iterations are distributed to the threads in chunks.
     * When a thread finishes its chunk, it requests a new one.
     * The chunk size remains constant.
     */
    omp_dynamic,

    /**
     * Standard OpenMP guided: iterations are distributed to the threads in chunks.
     * When a thread finishes its chunk, it requests a new one.
     * The chunk size starts large, and decreases over time towards the minimum set by the chunk size argument.
     */
    omp_guided,

    /**
     * Runtime: uses the scheduling kind set by the OMP_SCHEDULE environment variable.
     */
    omp_runtime,

    /**
     * Standard OpenMP static: iterations are distributed to the threads in chunks.
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
     * ExhaustiveSel: tries each algorithm once, then selects the fastest one and keeps it for a while.
     * Load imbalance re-triggers selection.
     */
    auto4omp_exhaustivesel,

    /**
     * Binary Search: tries the extreme ends of the algorithm list (ordered by load balancing overhead),
     * moves to the half of the list where the better algorithm is, and searches there.
     * Repeats until converging to the best scheduling technique.
     */
    auto4omp_binarySearch,

    /**
     * ExpertSel: uses runtime performance info and fuzzy logic with expert rules to select scheduling algorithms.
     */
    auto4omp_expertsel,

#ifdef AUTOPAS_USE_LB4OMP
    // LB4OMP's scheduling techniques [1, 2] (beware, technique names in the papers and Git README are outdated):
    /**
     * Profiling: uses dynamic,1 and tracks the execution times of loop iterations.
     * Some other techniques use the gathered information in future runs.
     * Pass a path for storing the data with KMP_PROFILE_DATA.
     */
    lb4omp_profiling,

    /**
     * Fixed Size Chunk: calculates an optimal chunk size and uses it for the entire simulation.
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_fsc,

    /**
     * Modified Fixed Size Chunk: calculates an optimal chunk size and uses it for the entire simulation.
     * Does not require profiling.
     */
    lb4omp_mfsc,

    /**
     * Tapering: similar to guided, but uses a more optimal model to decrease the chunk size.
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_tap,

    /**
     * Factoring: schedules the chunks in batches to the threads.
     * Each batch uses a specific chunk size calculated by its first thread.
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_fac,

    /**
     * Improved Factoring: like fac, but each thread calculates the chunk size of the batch independently,
     * instead of idly waiting for the batch's first thread to calculate it.
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_faca,

    /**
     * Bold:  similar to fac, but starts with larger chunks to reduce the scheduling overhead.
     * Beware: requires a prior profiling run of the simulation. Point KMP_PROFILE_DATA to the profiling data file.
     */
    lb4omp_bold,

    /**
     * Practical Factoring: like fac, but new batches pack half the remaining iterations.
     */
    lb4omp_fac2,

    /**
     * Weighted Factoring: like fac, but the sizes of chunks in a same batch vary
     * based on the processing weight of their host thread.
     * The processing weights of threads represent the specs of the computing node that runs them.
     */
    lb4omp_wf,

    /**
     * Adaptive Factoring: like fac, but uses run-time info instead of profiling data.
     */
    lb4omp_af,

    /**
     * Adaptive Weighted Factoring: like wf, but processing weights change dynamically based on thread performance.
     */
    lb4omp_awf,

    /**
     * Trapezoid Factoring Self Scheduling: like fac, but the chunk size decreases linearly over time.
     */
    lb4omp_tfss,

    /**
     * Fixed Increase Self Scheduling: periodically increases the chunk size by a constant bump.
     */
    lb4omp_fiss,

    /**
     * Variable Increase Self Scheduling: periodically increases the chunk size by a variable bump.
     * The bump starts as half the initial chunk size, and gets halved each period.
     */
    lb4omp_viss,

    /**
     * Random:  the chunk size varies randomly within a specific range.
     */
    lb4omp_rnd,

    // LB4OMP's scheduling techniques used by Auto4OMP (in addition to the standard scheduling kinds) [1, 2]:
    /**
     * Trapezoid Self Scheduling (from standard OpenMP): similar to guided,
     * but the chunk size decreases linearly.
     */
    lb4omp_trapezoidal,

    /**
     * Improved Practical Factoring: fac2 plus the optimizations from faca.
     */
    lb4omp_fac2a,

    /**
     * Static with Steal enabled (from standard OpenMP): like static,
     * but threads steal tasks from each other if they finish early.
     */
    lb4omp_static_steal,

    /**
     * Adaptive Weighted Factoring Variant B: like awf, but processing weights are updated at the end of each batch.
     */
    lb4omp_awf_b,

    /**
     * Adaptive Weighted Factoring Variant C: like awf, but  processing weights are updated at the end of each chunk.
     */
    lb4omp_awf_c,

    /**
     * Adaptive Weighted Factoring Variant D: like awf_c, but scheduling overhead influences processing weights.
     */
    lb4omp_awf_d,

    /**
     * Adaptive Weighted Factoring Variant E: like awf_b, but scheduling overhead influences processing weights.
     */
    lb4omp_awf_e,

    /**
     * Improved Adaptive Factoring: like af, but prior scheduling overhead influences chunk size.
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
   * Converts the old LB4OMP scheduling technique names to their corresponding new names.
   * The old names are the ones used in the LB4OMP and Auto4OMP papers and Git Readme.
   * The new names are the ones accepted in practice by the master branch LB4OMP.
   * If no conversion takes place, the given name will be returned in lower case.
   * @param name the name to convert
   * @return the converted lower case name
   */
   static std::string toNewName(const std::string &name) {
     std::string lowName = name;
     std::transform(lowName.begin(), lowName.end(), lowName.begin(), ::tolower);
     if (autopas::utils::StringUtils::contains(lowName, "ss")) return "dynamic";
     else if (autopas::utils::StringUtils::contains(lowName, "gss")) return "guided";
     else if (autopas::utils::StringUtils::contains(lowName, "mfac2")) return "fac2a";
     else if (autopas::utils::StringUtils::contains(lowName, "mfac")) return "faca";
     else if (autopas::utils::StringUtils::contains(lowName, "awf-b")) return "awf_b";
     else if (autopas::utils::StringUtils::contains(lowName, "awf-c")) return "awf_c";
     else if (autopas::utils::StringUtils::contains(lowName, "awf-d")) return "awf_d";
     else if (autopas::utils::StringUtils::contains(lowName, "awf-e")) return "awf_e";
     else if (autopas::utils::StringUtils::contains(lowName, "tss")) return "trapezoidal";
     else if (autopas::utils::StringUtils::contains(lowName, "maf")) return "af_a";
     return lowName;
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
