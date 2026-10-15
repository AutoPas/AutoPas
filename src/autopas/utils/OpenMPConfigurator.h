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
#include "autopas/utils/OpenMPChunkSize.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * This class provides configurable parameters for OpenMP.
 */
class OpenMPConfigurator {
 private:
  /**
   * The OpenMP scheduling kind, LB4OMP scheduling technique, or Auto4OMP selection method to use.
   */
  OpenMPKindOption _kind;

  /**
   * The OpenMP scheduling chunk size to use.
   */
  OpenMPChunkSize _chunkSize;

  /**
   * Number of iterations of each loop the schedule is used for.
   */
  std::vector<size_t> _loopCounts{};

  /**
   * Number of threads of the team that executes the loops.
   */
  size_t _numThreads{1};

 public:
  /**
   * OpenMP configurator constructor.
   * @param kind the OpenMP scheduling kind, LB4OMP scheduling technique, or Auto4OMP selection method to use
   * @param chunkSize the OpenMP scheduling chunk size to use
   * @param loopCounts Number of iterations of each loop the schedule is used for, e.g. the loops of a traversal (see
   * TraversalInterface::getRuntimeScheduledLoopCounts()). For collapsed loops, this is the number of iterations of the
   * collapsed loop. If this is empty e.g. if they are not known: the golden chunk size defaults to 1, which the logic
   * handler eliminates as an invalid configuration.
   * @param numThreads Number of threads.
   */
  [[maybe_unused]] explicit OpenMPConfigurator(OpenMPKindOption kind, OpenMPChunkSize chunkSize,
                                               const std::vector<size_t> &loopCounts, size_t numThreads);

  /**
   * Returns LB4OMP's golden (in papers called "expert") chunk size for a loop (see goldenChunkSize() in LB4OMP's
   * kmp_dispatch.cpp or refer to [2] (see end of file)).
   *
   * We make a minor modification to the formula, using a chunk size of 1 where the formula would give a chunk size of
   * 0.
   *
   * @param loopCount Number of iterations of the loop. For collapsed loops, this is the number of iterations of the
   * collapsed loop.
   * @param numThreads Number of threads of the team that executes the loop.
   * @return The golden/expert chunk size.
   */
  [[maybe_unused]] [[nodiscard]] static size_t getGoldenChunkSize(size_t loopCount, size_t numThreads);

  /**
   * AutoPas OpenMP configurator chunk size getter.
   * @return the current OpenMP chunk size
   */
  [[maybe_unused]] [[nodiscard]] OpenMPChunkSize getChunkSize() const;

  /**
   * OpenMP chunk size getter for setting OpenMP's scheduling runtime variables. Converts from our convenience/clarity
   * helpers to the actually taken chunk sizes:
   * - The Auto4OMP automated scheduling techniques are set with an auto schedule option with different chunk sizes.
   *   Here we convert from our helper scheduling kind options for them to the chunk size.
   * - Our special chunk sizes "balanced" and "golden". balanced becomes 0, which is what the OpenMP runtimes take for
   *   the default chunk size. golden gets calculated and returned.
   *
   * Throws for a fixed chunk size of 0 and for the not-applicable chunk size.
   *
   * @return the actual OpenMP chunk size, directly usable in OpenMP's schedule setter
   */
  [[maybe_unused]] [[nodiscard]] int getActualOMPChunkSize() const;

  /**
   * AutoPas OpenMP configurator chunk size setter.
   * @param chunkSize the new chunk size to use
   */
  [[maybe_unused]] void setChunkSize(OpenMPChunkSize chunkSize);

  /**
   * AutoPas OpenMP configurator scheduling kind getter.
   * @return the current OpenMP scheduling kind
   */
  [[maybe_unused]] [[nodiscard]] OpenMPKindOption getKind() const;

  /**
   * OpenMP scheduling kind getter for setting OpenMP's scheduling runtime variables.
   * Throws if the compiler and OpenMP runtime AutoPas is built with do not support the kind (see
   * OpenMPKindOption::isSupportedAtScheduleLevel()), and for the sequential kind, which is not an OpenMP kind.
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
   * Returns true if the OpenMP runtime would silently change the schedule kind to another "fallback" kind for all loops
   * of the configurator (typically due to chunk size and number of threads being too large), i.e. in cases of colored
   * traversals, if the fall back would occur for all colours. This is not problematic in the sense of producing errors,
   * but results in AutoPas running what it thinks are different configurations but in reality are actually duplicate
   * configurations. This function is therefore useful to eliminate such duplicate configurations from the search space.
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
   * @return True if the OpenMP runtime would fall back to another kind for every loop. False if there are no loops.
   */
  [[maybe_unused]] [[nodiscard]] bool fallsBackToOtherKindForAllLoops() const;

  /**
   * Returns true if a special chunk size (balanced/golden) would simply evaluate to a chunk size of 1 for the loops of
   * the configurator. This is the case if
   * - the golden chunk size evaluates to 1
   * - with static-balanced, there is no more than one iteration per thread.
   *
   * This function is useful for eliminating redundant duplicate configurations.
   *
   * @return True if the chunk size falls back to a chunk size of 1 for every loop. False for fixed and not-applicable
   * chunk sizes.
   */
  [[maybe_unused]] [[nodiscard]] bool fallsBackToChunkSize1() const;
};  // class OpenMPConfigurator

/**
 * Sets OpenMP's runtime schedule from a given OpenMP configurator.
 * schedule(runtime) will then use them for the traversal in the concerned calling thread.
 * @param ompConfig the OpenMP configurator
 */
inline void autopas_set_schedule(autopas::OpenMPConfigurator ompConfig) {
  autopas_set_schedule(ompConfig.getOMPKind(), ompConfig.getActualOMPChunkSize());
}  // void autopas_set_schedule
}  // namespace autopas

/*
 * Sources:
 * [1] https://www.computer.org/csdl/journal/td/2022/04/09524500/1wpqIcNI6YM
 * [2] https://ieeexplore.ieee.org/document/9825675
 */
