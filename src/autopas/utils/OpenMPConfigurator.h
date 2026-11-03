/**
 * @file OpenMPConfigurator.h
 * @author MehdiHachicha
 * @date 12.03.2024
 */

#pragma once

#include <cstddef>
#include <optional>
#include <set>
#include <utility>
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
   * Returns the schedule the OpenMP runtime actually runs if this is different from the actually configured one, for
   * all loops.
   *
   * This may occur for several reasons, such as the chunk size and number of threads being too large that another
   * schedule must be used, or when a special chunk size is used.
   *
   * This is not problematic in the sense of producing errors, but results in AutoPas running what it thinks is a
   * different configuration to what is actually run. This can lead to duplicate trials of the same actually run
   * configuration, under different names, and means that configurations outputted by AutoPas, e.g. the tuning loggers,
   * may be misleading.
   *
   * The following fallbacks are tested
   *
   * (1) With only one thread (which includes builds without OpenMP), every schedule runs the same as static, 1.
   *
   * (2) Fallbacks on LLVM-based libomp/LB4OMP:
   * - static_steal is only used if numThreads > 1 and ceil(loopCount / chunk) >= numThreads. Otherwise, LB4OMP and
   *   libomp < 10 use static, default chunk size, and libomp >= 10 uses dynamic, at the given chunk size.
   * - guided/auto are only used if numThreads > 1 and (2 * chunk + 1) * numThreads < loopCount. Otherwise, dynamic is
   *   used, or static if there is only one thread.
   *
   * GCC's libgomp does not fall back from guided. This is determined from OpenMPKindOption::runtimeScheduleLevel.
   *
   * (3) On every runtime, a special chunk size (balanced/golden) falls back to the fixed chunk size it evaluates to:
   * - golden evaluates to the fixed chunk size getActualOMPChunkSize() calculates. Not for Auto4OMP's selection
   *   methods, for which golden only marks the method.
   * - with static-balanced, there is no more than one iteration per thread, it is a chunk size of 1.
   *
   * @note This function will not work perfectly without being more invasive and high maintainence, but matches a large
   * number of scenarios.
   *
   * @return The schedule kind and chunk size the OpenMP runtime falls back to. Nothing if it runs the configured
   * schedule, or if there are no loop counts (unless there are no loop counts and golden chunk size is used, in which
   * case we warn and provide a fallback chunk size of 1).
   */
  [[maybe_unused]] [[nodiscard]] std::optional<std::pair<OpenMPKindOption, OpenMPChunkSize>> getFallbackSchedule()
      const;
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
