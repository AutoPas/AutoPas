/**
 * @file OpenMPConfigurator.cpp
 * @author MehdiHachicha
 * @date 13.03.2024
 */

#include "autopas/utils/OpenMPConfigurator.h"

namespace autopas {
/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-chunk-size <int>
 */
int openMPDefaultChunkSize = 0;

/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-kind <kind>
 */
OpenMPKindOption openMPDefaultKind = OpenMPKindOption::omp_runtime;

/**
 * OpenMP configurator default constructor.
 */
[[maybe_unused]] autopas::OpenMPConfigurator::OpenMPConfigurator() = default;

/**
 * OpenMP configurator constructor.
 * @param s chunk size used in OpenMP's loop scheduling
 */
[[maybe_unused]] autopas::OpenMPConfigurator::OpenMPConfigurator(OpenMPKindOption kind, int chunkSize) {
  setKind(kind);
  setChunkSize(chunkSize);
}

/**
 * OpenMP chunk size getter.
 * @return the current OpenMP chunk size
 */
[[maybe_unused]] [[nodiscard]] int autopas::OpenMPConfigurator::getChunkSize() const {
  return _kind;
}

/**
 * OpenMP chunk size getter for setting OpenMP's scheduling runtime variables.
 * @return the current OpenMP chunk size
 */
[[maybe_unused]] [[nodiscard]] int autopas::OpenMPConfigurator::getOMPChunkSize() const {
  switch (_kind) {
    case OpenMPKindOption::omp_auto:
      return 1;
    case OpenMPKindOption::auto4omp_randomsel:
      return 2;
    case OpenMPKindOption::auto4omp_exhaustivesel:
      return 3;
    case OpenMPKindOption::auto4omp_expertsel:
      return 4;
    default:
      return std::max(1, _chunkSize);
  }
}

/**
 * OpenMP chunk size setter.
 * @param s the new chunk size to use
 */
[[maybe_unused]] void autopas::OpenMPConfigurator::setChunkSize(int chunkSize) { _chunkSize = chunkSize; }

/**
 * OpenMP scheduling kind getter.
 * @return the current OpenMP scheduling kind, directly usable as an argument for OpenMP's schedule-clause
 */
[[maybe_unused]] [[nodiscard]] OpenMPKindOption autopas::OpenMPConfigurator::getKind() const { return _kind; }

/**
 * OpenMP standard scheduling kind getter.
 * @return the current OpenMP scheduling kind, directly usable as an argument for OpenMP's schedule-clause
 */
[[maybe_unused]] [[nodiscard]] omp_sched_t autopas::OpenMPConfigurator::getOMPKind() const {
  switch (_kind) {
    case OpenMPKindOption::omp_dynamic:
      return omp_sched_dynamic;
    case OpenMPKindOption::omp_guided:
      return omp_sched_guided;
    case OpenMPKindOption::omp_static:
      return omp_sched_static;
    default:
      return omp_sched_auto;
  }
}

/**
 * OpenMP chunk size setter.
 * @param k the new scheduling kind to use
 */
[[maybe_unused]] void autopas::OpenMPConfigurator::setKind(OpenMPKindOption kind) { _kind = kind; }

/**
 * Tells whether the scheduling chunk size should be overwritten.
 * @return whether the scheduling chunk size should be overwritten
 */
[[maybe_unused]] [[nodiscard]] bool autopas::OpenMPConfigurator::overrideChunkSize() const { return _kind >= 1; }
}  // namespace autopas