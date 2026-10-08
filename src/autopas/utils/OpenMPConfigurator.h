/**
 * @file OpenMPConfigurator.h
 * @author MehdiHachicha
 * @date 12.03.2024
 */

#pragma once

#include <cstddef>
#include <set>
#include <vector>

#include "autopas/options/OpenMPKindOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * OpenMP default chunk size.
 * md-flexible: set via command-line option --openmp-chunk-size int
 */
extern int openMPDefaultChunkSize;

/**
 * OpenMP default scheduling kind.
 * md-flexible: set via command-line option --openmp-kind kind
 */
extern OpenMPKindOption openMPDefaultKind;

/**
 * This class provides configurable parameters for OpenMP.
 */
class OpenMPConfigurator {
 private:
  /**
   * The OpenMP scheduling kind, LB4OMP scheduling technique, or Auto4OMP selection method to use.
   */
  OpenMPKindOption _kind = openMPDefaultKind;

  /**
   * The OpenMP scheduling chunk size to use.
   */
  int _chunkSize = openMPDefaultChunkSize;

 public:
  /**
   * OpenMP configurator default constructor.
   */
  [[maybe_unused]] OpenMPConfigurator();

  /**
   * OpenMP configurator constructor.
   * @param kind the OpenMP scheduling kind, LB4OMP scheduling technique, or Auto4OMP selection method to use
   * @param chunkSize the OpenMP scheduling chunk size to use
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
   * Throws if the compiler and OpenMP runtime AutoPas is built with do not support the kind (see
   * OpenMPKindOption::isSupportedAtScheduleLevel()).
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
   * Returns true if the OpenMP runtime would silently change a loop's schedule kind to another "fallback" kind
   * (typically due to chunk size and number of threads being too large). This is not problematic in the sense of
   * producing errors, but results in AutoPas running what it thinks are different configurations but in reality are
   * actually duplicate configurations.
   *
   * E.g. at a high chunk size and thread count and low loop length, static-steal falls back to e.g. dynamic scheduling
   * but AutoPas may have already trialled dynamic scheduling at that chunk size. This can result in a large number of
   * duplicate configurations being run, which AutoPas is not aware of. This function therefore manually determines if
   * a fallback would happen.
   *
   * Fallbacks on LLVM-based libomp/LB4OMP:
   * - static_steal is only used if numThreads > 1 and ceil(loopCount / chunk) >= numThreads. Otherwise, LB4OMP and
   *   libomp < 10 use static, default chunk size, and libomp >= 10 uses dynamic, at the given chunk size.
   * - guided/auto are only used if numThreads > 1 and (2 * chunk + 1) * numThreads < loopCount. Otherwise, dynamic is
   *   used, or static if there is only one thread.
   *
   * GCC's libgomp does not fall back from guided. This is determined from OpenMPKindOption::runtimeScheduleLevel.
   *
   * @note This function will not work perfectly without being more invasive and high maintainence, but should help
   * prune some duplicate configurations.
   *
   * @param loopCount Number of iterations of the loop. For collapsed loops, this is the number of iterations of the
   * collapsed loop.
   * @param numThreads Number of threads of the team that executes the loop.
   * @return True if the OpenMP runtime would fall back to another kind.
   */
  [[maybe_unused]] [[nodiscard]] bool fallsBackToOtherKind(size_t loopCount, int numThreads) const;

  /**
   * Tells whether the OpenMP runtime would fall back to another kind for all loops with the given loop counts (see
   * fallsBackToOtherKind()), i.e. in cases of coloured traversals, if the fall back would occur for all colours.
   * @param loopCounts Number of iterations of each loop.
   * @param numThreads Number of threads of the team that executes the loops.
   * @return True if the OpenMP runtime would fall back to another kind for every loop. False if there are no loops.
   */
  [[maybe_unused]] [[nodiscard]] bool fallsBackToOtherKindForAllLoops(const std::vector<size_t> &loopCounts,
                                                                      int numThreads) const;
};  // class OpenMPConfigurator

/**
 * Sets OpenMP's runtime schedule from a given OpenMP configurator.
 * schedule(runtime) will then use them for the traversal in the concerned calling thread.
 * @param ompConfig the OpenMP configurator
 */
inline void autopas_set_schedule(autopas::OpenMPConfigurator ompConfig) {
  autopas_set_schedule(ompConfig.getOMPKind(), ompConfig.getOMPChunkSize());
}  // void autopas_set_schedule
}  // namespace autopas

/*
 * Sources:
 * [1] https://www.computer.org/csdl/journal/td/2022/04/09524500/1wpqIcNI6YM
 * [2] https://ieeexplore.ieee.org/document/9825675
 */
