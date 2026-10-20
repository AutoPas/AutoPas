/**
 * @file OpenMPConfigurator.cpp
 * @author MehdiHachicha
 * @date 13.03.2024
 */

#include "autopas/utils/OpenMPConfigurator.h"

#include <algorithm>
#include <cmath>
#include <ranges>

#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/logging/Logger.h"

namespace autopas {
[[maybe_unused]] OpenMPConfigurator::OpenMPConfigurator(OpenMPKindOption kind, OpenMPChunkSize chunkSize,
                                                        const std::vector<size_t> &loopCounts, size_t numThreads)
    : _loopCounts(loopCounts), _numThreads(numThreads) {
  setKind(kind);
  setChunkSize(chunkSize);
}

[[maybe_unused]] [[nodiscard]] size_t OpenMPConfigurator::getGoldenChunkSize(size_t loopCount, size_t numThreads) {
  const auto numThreadsNonZero = std::max<size_t>(numThreads, 1);
  const auto iterationsPerThread = loopCount / numThreadsNonZero;
  // LB4OMP's formula gives 0 with fewer than two iterations per thread, and is undefined with fewer than one.
  // In these cases, we just use a chunk size of 1.
  if (iterationsPerThread < 2) {
    return 1;
  }
  // Golden ratio = 0.618: choose a chunk size in the "middle" between 1 and loopCount / (2 * numThreads).
  const auto exponent = std::floor(std::log2(iterationsPerThread) * 0.618);
  const auto numChunksPerThread = static_cast<size_t>(std::pow(2., exponent + 1.));
  return loopCount / (numChunksPerThread * numThreadsNonZero);
}

/**
 * AutoPas OpenMP configurator chunk size getter.
 * @return the current OpenMP chunk size
 */
[[maybe_unused]] [[nodiscard]] OpenMPChunkSize OpenMPConfigurator::getChunkSize() const { return _chunkSize; }

[[maybe_unused]] [[nodiscard]] int OpenMPConfigurator::getActualOMPChunkSize() const {
  // If an auto kind, return the appropriate chunk size.
  switch (_kind) {
    case OpenMPKindOption::omp_auto:
      return 1;
    case OpenMPKindOption::auto4omp_randomSel:
      return 2;
    case OpenMPKindOption::auto4omp_exhaustiveSel:
      return 3;
    case OpenMPKindOption::auto4omp_binarySearch:
      return 4;
    case OpenMPKindOption::auto4omp_expertSel:
      return 5;
    default:
      break;
  }

  switch (_chunkSize.getType()) {
    case OpenMPChunkSize::Type::fixed:
      if (_chunkSize.getSize() == 0) {
        utils::ExceptionHandler::exception(
            "OpenMPConfigurator::getActualOMPChunkSize(): A fixed chunk size must be at least 1. For a default chunk "
            "size for static, use OpenMPChunkSize::balanced; for other schedulers use 1.");
      }
      return static_cast<int>(_chunkSize.getSize());
    case OpenMPChunkSize::Type::balanced:
      // OpenMP runtimes use a chunk size of 0 for the default chunk size i.e. in the case of static scheduling, this is
      // balanced
      return 0;
    case OpenMPChunkSize::Type::golden: {
      // We use the smallest golden chunk size of the golden chunk size for all loops.
      size_t goldenChunkSize = 1;
      if (not _loopCounts.empty()) {
        goldenChunkSize = std::ranges::min(_loopCounts | std::views::transform([&](const auto loopCount) {
                                             return getGoldenChunkSize(loopCount, _numThreads);
                                           }));
      }
      return static_cast<int>(goldenChunkSize);
    }
    case OpenMPChunkSize::Type::notApplicable:
      utils::ExceptionHandler::exception(
          "OpenMPConfigurator::getActualOMPChunkSize(): The not-applicable chunk size is only used with sequential "
          "traversals, for which no schedule should be set.");
      return 0;
  }
  utils::ExceptionHandler::exception("OpenMPConfigurator::getActualOMPChunkSize(): Unknown type of chunk size.");
  return 0;
}

/**
 * AutoPas OpenMP configurator chunk size setter.
 * @param chunkSize the new chunk size to use
 */
[[maybe_unused]] void OpenMPConfigurator::setChunkSize(OpenMPChunkSize chunkSize) { _chunkSize = chunkSize; }

/**
 * AutoPas OpenMP configurator scheduling kind getter.
 * @return the current OpenMP scheduling kind
 */
[[maybe_unused]] [[nodiscard]] OpenMPKindOption OpenMPConfigurator::getKind() const { return _kind; }

/**
 * OpenMP scheduling kind getter for setting OpenMP's scheduling runtime variables.
 * @return the current OpenMP kind, directly usable in OpenMP's schedule setter
 */
[[maybe_unused]] [[nodiscard]] omp_sched_t OpenMPConfigurator::getOMPKind() const {
  if (not _kind.isSupportedAtScheduleLevel(OpenMPKindOption::runtimeScheduleLevel)) {
    utils::ExceptionHandler::exception(
        "OpenMPConfigurator::getOMPKind(): The OpenMP schedule kind {} is not supported by the compiler and OpenMP "
        "runtime AutoPas is built with. Configuration::hasCompatibleValues() should have rejected it.",
        _kind.to_string());
  }

  switch (_kind) {
    case OpenMPKindOption::omp_dynamic:
      return omp_sched_dynamic;
    case OpenMPKindOption::omp_guided:
      return omp_sched_guided;
    case OpenMPKindOption::omp_static:
      return omp_sched_static;
    // libomp's extensions. Declared by LB4OMP's omp.h, otherwise by WrapOpenMP.h.
    case OpenMPKindOption::llvm_trapezoidal:
      return omp_sched_trapezoidal;
    case OpenMPKindOption::lb4omp_static_steal:
      return omp_sched_static_steal;
#if AUTOPAS_OPENMP_SCHEDULE_LEVEL >= 3  // OpenMPKindOption::ScheduleLevel::clangLB4OMP
    // LB4OMP's scheduling techniques, only declared by LB4OMP's omp.h.
    case OpenMPKindOption::lb4omp_profiling:
      return omp_sched_profiling;
    case OpenMPKindOption::lb4omp_fsc:
      return omp_sched_fsc;
    case OpenMPKindOption::lb4omp_mfsc:
      return omp_sched_mfsc;
    case OpenMPKindOption::lb4omp_tap:
      return omp_sched_tap;
    case OpenMPKindOption::lb4omp_fac:
      return omp_sched_fac;
    case OpenMPKindOption::lb4omp_faca:
      return omp_sched_faca;
    case OpenMPKindOption::lb4omp_bold:
      return omp_sched_bold;
    case OpenMPKindOption::lb4omp_fac2:
      return omp_sched_fac2;
    case OpenMPKindOption::lb4omp_wf:
      return omp_sched_wf;
    case OpenMPKindOption::lb4omp_af:
      return omp_sched_af;
    case OpenMPKindOption::lb4omp_awf:
      return omp_sched_awf;
    case OpenMPKindOption::lb4omp_tfss:
      return omp_sched_tfss;
    case OpenMPKindOption::lb4omp_fiss:
      return omp_sched_fiss;
    case OpenMPKindOption::lb4omp_fac2a:
      return omp_sched_fac2a;
    case OpenMPKindOption::lb4omp_awf_b:
      return omp_sched_awf_b;
    case OpenMPKindOption::lb4omp_awf_c:
      return omp_sched_awf_c;
    case OpenMPKindOption::lb4omp_awf_d:
      return omp_sched_awf_d;
    case OpenMPKindOption::lb4omp_awf_e:
      return omp_sched_awf_e;
    case OpenMPKindOption::lb4omp_af_a:
      return omp_sched_af_a;
#endif
    case OpenMPKindOption::sequential:
      utils::ExceptionHandler::exception(
          "OpenMPConfigurator::getOMPKind(): The sequential kind is not an OpenMP schedule kind. It is only used with "
          "sequential traversals, for which no schedule should be set.");
      return omp_sched_static;
    default:
      return omp_sched_auto;  // Standard auto and Auto4OMP's selection methods.
  }
}

/**
 * AutoPas OpenMP configurator scheduling kind setter.
 * @param kind the new scheduling kind to use
 */
[[maybe_unused]] void OpenMPConfigurator::setKind(OpenMPKindOption kind) { _kind = kind; }

[[maybe_unused]] [[nodiscard]] std::optional<std::pair<OpenMPKindOption, OpenMPChunkSize>>
OpenMPConfigurator::getFallbackSchedule() const {
  // Sequential traversals have no schedule.
  if (_kind == OpenMPKindOption::sequential) {
    return std::nullopt;
  }

  // (1) With one thread, every schedule runs the same as static, 1.
  if (_numThreads <= 1) {
    if (_kind == OpenMPKindOption::omp_static and _chunkSize == OpenMPChunkSize(1)) {
      return std::nullopt;
    }
    return {{OpenMPKindOption::omp_static, 1}};
  }

  // Without loop counts, nothing is known about the fallbacks below.
  if (_loopCounts.empty()) {
    if (_chunkSize == OpenMPChunkSize::golden and not OpenMPKindOption::getAuto4OMPSelectionMethods().contains(_kind)) {
      AutoPasLog(WARN,
                 "OpenMPConfigurator: The golden chunk size is used with a traversal that provides no loop counts, so "
                 "it falls back to a chunk size of 1. Either modify this traversal to provide loop counts, or remove "
                 "the golden chunk size from your search space.");
      return {{_kind, 1}};
    }
    return std::nullopt;
  }

  // (2) Fallbacks of the kind. GCC's libgomp has none. For LLVM's libomp, The kinds fall back for loops which are too
  // short, so all loops fall back if the longest one does.
  if (OpenMPKindOption::runtimeScheduleLevel != OpenMPKindOption::ScheduleLevel::gccLibgomp) {
    const auto longestLoopCount = std::ranges::max(_loopCounts);
    const auto chunkSize = static_cast<size_t>(getActualOMPChunkSize());

    switch (_kind) {
      case OpenMPKindOption::lb4omp_static_steal: {
        // static_steal needs at least one chunk per thread.
        const auto numChunks = longestLoopCount / chunkSize + (longestLoopCount % chunkSize == 0 ? 0 : 1);
        if (numChunks < _numThreads) {
          if (OpenMPKindOption::runtimeScheduleLevel == OpenMPKindOption::ScheduleLevel::clangLB4OMP) {
            return {{OpenMPKindOption::omp_static, OpenMPChunkSize::balanced}};
          }
          return {{OpenMPKindOption::omp_dynamic, _chunkSize}};
        }
        break;
      }
      case OpenMPKindOption::omp_guided:
      case OpenMPKindOption::omp_auto:
        // guided and auto need more than (2 * chunkSize + 1) iterations per thread. auto is a guided variant with a
        // chunk size of 1, regardless of the configuration chunk size.
        if ((2 * chunkSize + 1) * _numThreads >= longestLoopCount) {
          return {
              {OpenMPKindOption::omp_dynamic, _kind == OpenMPKindOption::omp_auto ? OpenMPChunkSize(1) : _chunkSize}};
        }
        break;
      default:
        break;
    }
  }

  // (3) Fallbacks of the special chunk sizes, on every runtime.
  switch (_chunkSize.getType()) {
    case OpenMPChunkSize::Type::golden:
      if (not OpenMPKindOption::getAuto4OMPSelectionMethods().contains(_kind)) {
        return {{_kind, static_cast<size_t>(getActualOMPChunkSize())}};
      }
      break;
    case OpenMPChunkSize::Type::balanced:
      // Static gives each thread one chunk of about loopCount / numThreads iterations. With at most one iteration per
      // thread, this is the same as a chunk size of 1.
      if (std::ranges::all_of(_loopCounts, [&](const auto loopCount) { return loopCount <= _numThreads; })) {
        return {{_kind, 1}};
      }
      break;
    case OpenMPChunkSize::Type::fixed:
    case OpenMPChunkSize::Type::notApplicable:
      break;
  }
  return std::nullopt;
}
}  // namespace autopas