/**
 * @file OpenMPConfigurator.h
 * @author MehdiHachicha
 * @date 12.03.2024
 */

#pragma once

#include <cstddef>
#include <set>

#include "autopas/options/Option.h"
#include "WrapOpenMP.h"

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
  };

  /**
   * Constructor.
   */
  OpenMPKindOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr OpenMPKindOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<OpenMPKindOption> getDiscouragedOptions() {
    return {};
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
        {OpenMPKindOption::omp_runtime, "runtime"},
        {OpenMPKindOption::omp_static, "static"},

        // Auto4OMP's automated selection methods:
        {OpenMPKindOption::auto4omp_randomsel, "randomSel"},
        {OpenMPKindOption::auto4omp_exhaustivesel, "exhaustiveSel"},
        {OpenMPKindOption::auto4omp_expertsel, "expertSel"},
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
 * md-flexible: set via command-line option --openmp-kind <int>
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
         * OpenMP schedule setter.
         * @param k the new scheduling kind
         * @param s the new chunk size
         */
        [[maybe_unused]] void setSchedule();

        /**
         * OpenMP schedule setter.
         * @param k the new scheduling kind
         * @param s the new chunk size
         */
        [[maybe_unused]] void setSchedule(OpenMPKindOption kind, int chunkSize);

        /**
         * Tells whether the scheduling chunk size should be overwritten.
         * @return whether the scheduling chunk size should be overwritten
         */
         [[maybe_unused]] [[nodiscard]] bool overrideChunkSize() const;
};
}