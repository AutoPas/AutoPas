/**
 * @file OpenMPConfigurator.cpp
 * @author MehdiHachicha
 * @date 13.03.2024
 */

#include "autopas/utils/OpenMPConfigurator.h"

namespace autopas {
/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-chunk-size <long>
 */
int openMPDefaultChunkSize = 0;

/**
 * OpenMP default chunk size for manual testing.
 * md-flexible: set via command-line option --openmp-chunk-size <long>
 */
OpenMPKindOption openMPDefaultKind = OpenMPKindOption::auto4omp_expertsel;

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
 * Sets OpenMP's scheduling runtime variables.
 */
[[maybe_unused]] inline void autopas::OpenMPConfigurator::setSchedule() const {
  autopas_set_schedule(getOMPKind(), getOMPChunkSize());
}

/**
 * Sets OpenMP's scheduling runtime variables.
 * @param k the new scheduling kind
 * @param s the new chunk size
 */
[[maybe_unused]] inline void autopas::OpenMPConfigurator::setSchedule(OpenMPKindOption kind, int chunkSize) {
  setKind(kind);
  setChunkSize(chunkSize);
  setSchedule();
}

/**
 * Gets OpenMP's scheduling runtime variables as a string.
 * @return the scheduling kind and chunk size used by schedule(runtime) as a string
 */
[[maybe_unused]] std::string autopas::OpenMPConfigurator::getRuntimeSchedule() {
  omp_sched_t kind;
  int chunkSize;
  autopas_get_schedule(&kind, &chunkSize);
  std::string schedule = "Schedule = ";
  switch (kind) {
    case omp_sched_static:
      schedule += "Static";
      break;
    case omp_sched_dynamic:
      schedule += "Dynamic";
      break;
    case omp_sched_guided:
      schedule += "Guided";
      break;
    case omp_sched_auto:
      schedule += "Auto";
      break;
    case omp_sched_monotonic:
      schedule += "Monotonic";
  }
  schedule += "," + std::to_string(chunkSize);
  return schedule;
}

/**
 * Tells whether the scheduling chunk size should be overwritten.
 * @return whether the scheduling chunk size should be overwritten
 */
[[maybe_unused]] [[nodiscard]] bool autopas::OpenMPConfigurator::overrideChunkSize() const { return _kind >= 1; }
}  // namespace autopas