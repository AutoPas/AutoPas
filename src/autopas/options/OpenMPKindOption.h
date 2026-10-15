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
   *
   * Which of these the OpenMP runtime supports is described by ScheduleLevel.
   *
   * Some kinds are deliberately not options:
   * - OpenMP's runtime kind, which takes the kind from the OMP_SCHEDULE environment variable: AutoPas always sets the
   *   schedule itself, so it makes no sense as a tunable option.
   * - LB4OMP's Variable Increase Self Scheduling (viss) and Random (rnd): LB4OMP's omp_get_schedule() does not handle
   *   them and aborts with KMP_FATAL. These could be fixed in LB4OMP, but doubt they will be effective in AutoPas.
   */
  enum Value {
    // Standard OpenMP's scheduling kinds:
    /**
     * OpenMP static scheduling.
     *
     * Iterations are distributed to the threads at the beginning of the loop. Where a chunk size is provided, each
     * chunk is distributed in a round-robin fashion.
     *
     * In AutoPas, static scheduling can be paired with a balanced chunk size, resulting in one chunk per thread that
     * evenly divides the loop between them (the default static behavior, where no chunk size is provided).
     *
     * Pros: No scheduling overhead, and predictable memory access.
     * Cons: No ability to dynamically load balance, if some chunks take more time.
     *
     * Smaller chunks = better load balance (lower thread idle time). Large chunks/balanced = better memory access.
     *
     * Both static with a balanced chunk size and with a chunk size of 1 are good choices for AutoPas. Don't
     * underestimate the ability of static with a chunk size of 1 to load balance.
     */
    omp_static,

    /**
     * OpenMP dynamic scheduling.
     *
     * Chunks are distributed to the threads once the thread is done with their previous chunk.
     *
     * Pros: No predetermined distribution of loops allows better ability to adapt to load imbalance => better load
     * balance, lower thread idle times than static.
     * Cons: Unpredictable memory access, some scheduling overhead (typically minor with AutoPas).
     *
     * Larger chunk sizes reduce the impact of the pros but also the impact of the cons.
     *
     * dynamic with a chunk size of 1 is a good choice for AutoPas in heterogeneous scenarios, although static_steal is
     * typically better in such scenarios, but is not available on GCC/libgomp and older Clang/libomp.
     */
    omp_dynamic,

    /**
     * OpenMP guided scheduling.
     *
     * Like dynamic, but begin with larger chunk sizes and reduce them to a minimum chunk size (that provided by us).
     *
     * Pros: combines the pros of dynamic with a small chunk size and the pros of dynamic with a large chunk size.
     * Cons: combines the cons of dynamic with a small chunk size and the cons of dynamic with a large chunk size.
     *
     * Falls back to dynamic for short loops (see OpenMPConfigurator::fallsBackToOtherKind()).
     *
     * guided is rarely found to be that good a choice for AutoPas, we rarely have large enough loops.
     */
    omp_guided,

    /**
     * OpenMP auto scheduling.
     *
     * In practice, never really implemented as an actual automated schedule selection.
     *
     * For GCC/libgomp, it is implemented just as static with the default chunk size (i.e. static with a balanced chunk
     * size), so we just skip this when trialling configurations. For Clang/libomp/LB4OMP, it is implemented as an analytical variant of guided.
     *
     * See guided for pros and cons.
     *
     * Similarly, falls back to dynamic for short loops (see OpenMPConfigurator::fallsBackToOtherKind()).
     *
     * This can compare well to standard guided, but is similarly rarely found to be a good choice for AutoPas.
     */
    omp_auto,

    // Sequential
    /**
     * Sequential: not an OpenMP scheduling kind, but used with sequential traversals. It is never set by the user, but
     * implicitly included in the search space.
     */
    sequential,

    // LLVM's Intel extensions
    // Intel added the following two scheduling kinds, not accessible by direct #pragma omp for schedule(...), but can
    // be set in the setter for schedule(runtime) loops.
    /**
     * LLVM libomp's trapezoidal scheduling.
     *
     * Similar to guided, but with a linear reduction of chunk size. This is similarly rarely found to be a good choice
     * for AutoPas.
     */
    llvm_trapezoidal,

    /**
     * LLVM libomp's static steal scheduling.
     *
     * Begin with a static balanced scheduling. When a thread finishes, it can "steal" work from other threads.
     *
     * Pros: The memory access benefits of static, balanced scheduling, with the ability to dynamically load balance.
     * Cons: The load balancing has typically more overhead than for dynamic/guided/similar schedules. For perfectly
     * balanced distributions, there is some overhead which static, balanced does not have.
     *
     * There is a bug in LLVM's libomp OpenMP runtime library, which means this schedule can't be set yet. A fix should
     * come in LLVM 24. Otherwise, this schedule can be used when compiling LB4OMP as the OpenMP runtime library.
     *
     * For short loops, static_steal falls back to either static, balanced or dynamic (depends on version), so we skip
     * this schedule in tuning when that would occur (see OpenMPConfigurator::fallsBackToOtherKind()).
     *
     * Unlike static, a non-default chunk size does not result in round-robin loop assignment -> the assignment is
     * always balanced. The chunk size affects the granularity of the loops for the initial distribution and for
     * stealing.
     *
     * AutoPas performs excellently with static_steal and a chunk size of 1. Typically, the only cases where it is
     * outperformed are homogeneous ones (static, balanced is better) or for very sparse ones (where dynamic is better).
     */
    lb4omp_static_steal,

    // AutoPas can also be compiled with LB4OMP, which includes a fix for static_steal, an extended range of scheduling
    // techniques, and proper auto schedule selection methods (Auto4OMP). See [1, 2].

    // LB4OMP's scheduling techniques used by Auto4OMP:

    /**
     * LB4OMP's improved practical factoring scheduling (mFAC2 in [1, 2]).
     *
     * Chunks are handed out dynamically in batches. Each batch is split evenly per thread. The first batch covers
     * half of the loop, each following batch covers half of the previous, down to a minimum chunk size (that's provided
     * by us). Each thread calculates the chunk sizes itself, without synchronizing with the other threads.
     *
     * Performs similarly to guided, and only extremely rarely is optimal in AutoPas
     */
    lb4omp_fac2a,

    /**
     * LB4OMP's adaptive weighted factoring variant B scheduling (AWF-B in [1, 2]).
     *
     * Like fac2a, but each chunk size is multiplied by a weight of the thread that executes it. Whenever a thread
     * requests a new chunk, it updates its weight from the measured execution times of its previous chunks, relative to
     * those of the other threads, so that faster threads get larger chunks. Handles really more hardware level imbalance
     * than simulation scenario imbalance.
     *
     * Performs worse than fac2a for AutoPas. We don't really have enough loops for this to pay off.
     */
    lb4omp_awf_b,

    /**
     * LB4OMP's adaptive weighted factoring variant C scheduling (AWF-C in [1, 2]).
     *
     * Like awf_b, but the chunk size before weighting is recalculated for every chunk as the remaining iterations
     * divided by 2P, instead of once per batch.
     *
     * Performs worse than fac2a for AutoPas. We don't really have enough loops for this to pay off.
     */
    lb4omp_awf_c,

    /**
     * LB4OMP's adaptive weighted factoring variant D scheduling (AWF-D in [1, 2]).
     *
     * Like awf_b, but the measured times include the scheduling overhead, so that it also influences the weights.
     *
     * Performs worse than fac2a for AutoPas. We don't really have enough loops for this to pay off.
     */
    lb4omp_awf_d,

    /**
     * LB4OMP's adaptive weighted factoring variant E scheduling (AWF-E in [1, 2]).
     *
     * Like awf_c, but the measured times include the scheduling overhead, so that it also influences the weights.
     *
     * Performs worse than fac2a for AutoPas. We don't really have enough loops for this to pay off.
     */
    lb4omp_awf_e,

    /**
     * LB4OMP's improved adaptive factoring scheduling (mAF in [1, 2]).
     *
     * Each thread measures the mean and standard deviation of its iteration execution times at runtime. Each chunk size
     * is then calculated from these statistics of all threads and the number of remaining iterations, following the
     * factoring model, down to a minimum chunk size (that provided by us). To measure the statistics, every chunk is
     * executed as ten sub-chunks, and the first chunk of each thread has 100 iterations. The measured times include the
     * scheduling overhead (unlike af).
     *
     * Performs worse than fac2a for AutoPas. We don't really have enough loops for this to pay off.
     */
    lb4omp_af_a,

    // Auto4OMP's automated selection methods:
    /**
     * Auto4OMP's RandomSel automated schedule selector.
     *
     * May reselect a random scheduling algorithm at the end of each loop. The probability of reselection increases with
     * greater load imbalance.
     */
    auto4omp_randomSel,

    /**
     * Auto4OMP's ExhaustiveSel automated schedule selector.
     *
     * Tries each algorithm once, then selects the fastest one and keeps it for a while. Similar to AutoPas's Full
     * Search. Load imbalance re-triggers selection.
     */
    auto4omp_exhaustiveSel,

    /**
     * Auto4OMP's Binary Search automated schedule selector.
     *
     * Tries the extreme ends of the Auto4OMP scheduling algorithm list (ordered by load balancing overhead, see [2]),
     * moves to the half of the list where the better algorithm is, and searches there. Repeats until converging to the
     * best scheduling technique. Load imbalance re-triggers selection.
     */
    auto4omp_binarySearch,

    /**
     * Auto4OMP's ExpertSel automated schedule selector.
     *
     * Uses runtime performance info and fuzzy logic with expert rules to select scheduling algorithms and change them
     * as the performance info changes.
     */
    auto4omp_expertSel,

    // LB4OMP's remaining scheduling techniques [1, 2] (beware, technique names in the papers and Git README are outdated):
    // We do not recommend using these as they are untested by either tests or in our own experiments.
    /**
     * LB4OMP's profiling "scheduling".
     *
     * Schedules like dynamic with a chunk size of 1, and measures the mean and standard deviation of the iteration
     * execution times of each loop. These are stored in the file given by KMP_PROFILE_DATA, from which fsc, tap, fac,
     * faca, and bold read them in later runs.
     */
    lb4omp_profiling,

    /**
     * LB4OMP's fixed size chunking scheduling (FSC in [1]).
     *
     * Like dynamic, but with a chunk size calculated to balance load imbalance against scheduling overhead, from the
     * standard deviation of the iteration execution times and the scheduling overhead (given by KMP_OVERHEAD).
     *
     * Requires a prior profiling run (see lb4omp_profiling).
     */
    lb4omp_fsc,

    /**
     * LB4OMP's modified fixed size chunking scheduling (mFSC in [1]).
     *
     * Like fsc, but the chunk size is calculated only from the number of iterations and threads, so no profiling run is
     * needed.
     */
    lb4omp_mfsc,

    /**
     * LB4OMP's tapering scheduling (TAP in [1]).
     *
     * Like guided, but the chunk sizes decrease according to a model that also takes the mean and standard deviation of
     * the iteration execution times into account (scaled by KMP_ALPHA).
     *
     * Requires a prior profiling run (see lb4omp_profiling).
     */
    lb4omp_tap,

    /**
     * LB4OMP's factoring scheduling (FAC in [1]).
     *
     * Like fac2a, chunks are handed out in batches of one chunk per thread. The chunk size of each batch is calculated
     * from the number of remaining iterations and the mean and standard deviation of the iteration execution times, by
     * the thread that starts the batch, in a critical section.
     *
     * Requires a prior profiling run (see lb4omp_profiling).
     */
    lb4omp_fac,

    /**
     * LB4OMP's improved factoring scheduling (mFAC in [1]).
     *
     * Like fac, but each thread calculates the chunk sizes itself, instead of the thread that starts a batch doing so
     * in a critical section.
     *
     * Requires a prior profiling run (see lb4omp_profiling).
     */
    lb4omp_faca,

    /**
     * LB4OMP's BOLD scheduling.
     *
     * Like fac, but starts with larger chunks to reduce the scheduling overhead, and adapts the chunk sizes to the
     * progress of the threads. Also uses the scheduling overhead (given by KMP_OVERHEAD).
     *
     * Requires a prior profiling run (see lb4omp_profiling).
     */
    lb4omp_bold,

    /**
     * LB4OMP's practical factoring scheduling (FAC2 in [1]).
     *
     * The same chunk sizes as fac2a, calculated with a closed formula instead of by halving the chunk size of the
     * previous batch.
     */
    lb4omp_fac2,

    /**
     * LB4OMP's weighted factoring scheduling (WF2 in [1]).
     *
     * Like fac2a, but each chunk size is multiplied by a fixed weight of the thread that executes it, given by
     * KMP_WEIGHTS, e.g. for heterogeneous cores. Without KMP_WEIGHTS, all weights are 1, which is identical to fac2a.
     */
    lb4omp_wf,

    /**
     * LB4OMP's adaptive factoring scheduling (AF in [1]).
     *
     * Like af_a, but the measured times exclude the scheduling overhead, and the number of sub-chunks is given by
     * KMP_DIVIDER.
     */
    lb4omp_af,

    /**
     * LB4OMP's adaptive weighted factoring scheduling for time-stepping applications (AWF in [1]).
     *
     * Like wf, but the weights are updated at the end of every execution of the loop (i.e. every time step), from the
     * measured execution times of the threads. Unlike awf_b to awf_e, the weights stay fixed within one execution of the
     * loop.
     */
    lb4omp_awf,

    /**
     * LB4OMP's trapezoid factoring self scheduling (TFSS in [1]).
     *
     * Like trapezoidal, the chunk sizes decrease linearly, but chunks are handed out in batches of one chunk per thread.
     * The chunk size of a batch is the average of the next P chunk sizes of trapezoidal, for P threads.
     */
    lb4omp_tfss,

    /**
     * LB4OMP's fixed increase self scheduling (FISS in [1]).
     *
     * Chunks are handed out in batches of one chunk per thread. The chunk size starts small and increases by a fixed
     * amount from batch to batch. The number of increases is given by KMP_FISS_STAGES, by default P - 1 for P threads.
     */
    lb4omp_fiss,
  };

  /**
   * Levels of OpenMP schedule kind support, depending on the compiler and OpenMP runtime AutoPas is built with.
   * Each level supports the kinds of all lower levels, so levels can be compared with inequalities.
   * The level AutoPas is built with is runtimeScheduleLevel, which CMake sets (see autopas_OpenMP.cmake).
   */
  enum class ScheduleLevel : int {
    /**
     * GCC + libgomp: static, dynamic, guided. auto is implemented as static with default chunk size, so it is
     * considered unsupported to avoid redundantly duplicating static.
     */
    gccLibgomp = 0,
    /**
     * Clang + libomp: adds auto (implemented as a variant of guided) and trapezoidal. Whilst static_steal is
     * implemented, the setter does not work with it, so this is excluded.
     */
    clangLibomp = 1,
    /**
     * Clang + a libomp with a working static_steal: adds static_steal. This is a placeholder for now, until a fix
     * gets added to a later version of LLVM.
     */
    clangFixedLibomp = 2,
    /**
     * Clang + LB4OMP (Auto4OMP): adds all remaining Auto4OMP selection methods and LB4OMP scheduling techniques.
     */
    clangLB4OMP = 3,
  };

  static_assert(AUTOPAS_OPENMP_SCHEDULE_LEVEL >= static_cast<int>(ScheduleLevel::gccLibgomp) and
                    AUTOPAS_OPENMP_SCHEDULE_LEVEL <= static_cast<int>(ScheduleLevel::clangLB4OMP),
                "AUTOPAS_OPENMP_SCHEDULE_LEVEL is not a valid OpenMPKindOption::ScheduleLevel.");

  /**
   * The schedule level of the compiler and OpenMP runtime AutoPas is built with.
   */
  static constexpr ScheduleLevel runtimeScheduleLevel = static_cast<ScheduleLevel>(AUTOPAS_OPENMP_SCHEDULE_LEVEL);

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
   * Whether the given schedule level supports this kind, i.e. whether a compiler and OpenMP runtime with this level can
   * set it.
   * @param level The schedule level. For the level AutoPas is built with, use runtimeScheduleLevel.
   * @return True if the kind is supported at the given level.
   */
  [[nodiscard]] bool isSupportedAtScheduleLevel(ScheduleLevel level) const {
    switch (_value) {
      case omp_static:
      case omp_dynamic:
      case omp_guided:
      case sequential:
        return level >= ScheduleLevel::gccLibgomp;
      case omp_auto:
      case llvm_trapezoidal:
        return level >= ScheduleLevel::clangLibomp;
      case lb4omp_static_steal:
        return level >= ScheduleLevel::clangFixedLibomp;
      case auto4omp_randomSel:
      case auto4omp_exhaustiveSel:
      case auto4omp_binarySearch:
      case auto4omp_expertSel:
      case lb4omp_profiling:
      case lb4omp_fsc:
      case lb4omp_mfsc:
      case lb4omp_tap:
      case lb4omp_fac:
      case lb4omp_faca:
      case lb4omp_bold:
      case lb4omp_fac2:
      case lb4omp_wf:
      case lb4omp_af:
      case lb4omp_awf:
      case lb4omp_tfss:
      case lb4omp_fiss:
      case lb4omp_fac2a:
      case lb4omp_awf_b:
      case lb4omp_awf_c:
      case lb4omp_awf_d:
      case lb4omp_awf_e:
      case lb4omp_af_a:
        return level >= ScheduleLevel::clangLB4OMP;
      default:
        utils::ExceptionHandler::exception("OpenMPKindOption::isSupportedAtScheduleLevel(): Unknown option {}.",
                                           static_cast<int>(_value));
    }
    return false;
  }

  /**
   * Set of options that are very unlikely to be interesting.
   *
   * Among them:
   * - Auto4OMP's selection methods pick a scheduling technique themselves, which overlaps with what AutoPas' tuning
   *   does, and are currently broken.
   * - Several of LB4OMP's techniques need data that LB4OMP only has after a profiling run (see KMP_PROFILE_DATA), and
   *   abort inside the runtime otherwise.
   *
   * Also includes sequential, which is never set by the user, as it is always included implicitly.
   *
   * @return
   */
  static std::set<OpenMPKindOption> getDiscouragedOptions() {
    return {
        OpenMPKindOption::sequential,
        OpenMPKindOption::auto4omp_randomSel,
        OpenMPKindOption::auto4omp_exhaustiveSel,
        OpenMPKindOption::auto4omp_binarySearch,
        OpenMPKindOption::auto4omp_expertSel,
        OpenMPKindOption::lb4omp_profiling,
        OpenMPKindOption::lb4omp_fsc,
        OpenMPKindOption::lb4omp_mfsc,
        OpenMPKindOption::lb4omp_tap,
        OpenMPKindOption::lb4omp_fac,
        OpenMPKindOption::lb4omp_faca,
        OpenMPKindOption::lb4omp_bold,
        OpenMPKindOption::lb4omp_fac2,
        OpenMPKindOption::lb4omp_wf,
        OpenMPKindOption::lb4omp_af,
        OpenMPKindOption::lb4omp_awf,
        OpenMPKindOption::lb4omp_tfss,
        OpenMPKindOption::lb4omp_fiss,
        OpenMPKindOption::lb4omp_awf_b,
        OpenMPKindOption::lb4omp_awf_c,
        OpenMPKindOption::lb4omp_awf_d,
        OpenMPKindOption::lb4omp_awf_e,
        OpenMPKindOption::lb4omp_af_a,
    };
  }

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
        {OpenMPKindOption::omp_static, "static"},

        // Auto4OMP's automated selection methods:
        {OpenMPKindOption::auto4omp_randomSel, "randomSel"},
        {OpenMPKindOption::auto4omp_exhaustiveSel, "exhaustiveSel"},
        {OpenMPKindOption::auto4omp_binarySearch, "binarySearch"},
        {OpenMPKindOption::auto4omp_expertSel, "expertSel"},

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

        // LB4OMP's scheduling techniques used by Auto4OMP (in addition to the standard scheduling kinds):
        {OpenMPKindOption::llvm_trapezoidal, "trapezoidal"},    // Trapezoid Self Scheduling
        {OpenMPKindOption::lb4omp_fac2a, "fac2a"},                // Improved Practical Factoring
        {OpenMPKindOption::lb4omp_static_steal, "static_steal"},  // Static with Steal enabled
        {OpenMPKindOption::lb4omp_awf_b, "awf_b"},                // Adaptive Weighted Factoring Variant B
        {OpenMPKindOption::lb4omp_awf_c, "awf_c"},                // Adaptive Weighted Factoring Variant C
        {OpenMPKindOption::lb4omp_awf_d, "awf_d"},                // Adaptive Weighted Factoring Variant D
        {OpenMPKindOption::lb4omp_awf_e, "awf_e"},                // Adaptive Weighted Factoring Variant E
        {OpenMPKindOption::lb4omp_af_a, "af_a"},                  // Improved Adaptive Factoring

        // AutoPas' kind for sequential traversals:
        {OpenMPKindOption::sequential, "sequential"},
    };
  };

  /**
   * Set of all options, excluding sequential.
   * @return
   */
  static std::set<OpenMPKindOption> getAllParallelKinds() {
    auto s = getAllOptions();
    s.erase(sequential);
    return s;
  }

  /**
   * Set of Auto4OMP's selection methods.
   * @return
   */
  static std::set<OpenMPKindOption> getAuto4OMPSelectionMethods() {
    return {auto4omp_randomSel, auto4omp_exhaustiveSel, auto4omp_binarySearch, auto4omp_expertSel};
  }

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
    if (autopas::utils::StringUtils::contains(lowName, "ss"))
      return "dynamic";
    else if (autopas::utils::StringUtils::contains(lowName, "gss"))
      return "guided";
    else if (autopas::utils::StringUtils::contains(lowName, "mfac2"))
      return "fac2a";
    else if (autopas::utils::StringUtils::contains(lowName, "mfac"))
      return "faca";
    else if (autopas::utils::StringUtils::contains(lowName, "awf-b"))
      return "awf_b";
    else if (autopas::utils::StringUtils::contains(lowName, "awf-c"))
      return "awf_c";
    else if (autopas::utils::StringUtils::contains(lowName, "awf-d"))
      return "awf_d";
    else if (autopas::utils::StringUtils::contains(lowName, "awf-e"))
      return "awf_e";
    else if (autopas::utils::StringUtils::contains(lowName, "tss"))
      return "trapezoidal";
    else if (autopas::utils::StringUtils::contains(lowName, "maf"))
      return "af_a";
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
