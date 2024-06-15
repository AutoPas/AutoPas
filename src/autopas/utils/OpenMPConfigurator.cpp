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
 * AutoPas OpenMP configurator chunk size getter.
 * @return the current OpenMP chunk size
 */
[[maybe_unused]] [[nodiscard]] int autopas::OpenMPConfigurator::getChunkSize() const {
  return _chunkSize;
}

/**
 * OpenMP chunk size getter for setting OpenMP's scheduling runtime variables.
 * @return the current OpenMP chunk size, directly usable in OpenMP's schedule setter
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
 * AutoPas OpenMP configurator chunk size setter.
 * @param chunkSize the new chunk size to use
 */
[[maybe_unused]] void autopas::OpenMPConfigurator::setChunkSize(int chunkSize) { _chunkSize = chunkSize; }

/**
 * AutoPas OpenMP configurator scheduling kind getter.
 * @return the current OpenMP scheduling kind
 */
[[maybe_unused]] [[nodiscard]] OpenMPKindOption autopas::OpenMPConfigurator::getKind() const { return _kind; }

/**
 * OpenMP scheduling kind getter for setting OpenMP's scheduling runtime variables.
 * @return the current OpenMP kind, directly usable in OpenMP's schedule setter
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
 * AutoPas OpenMP configurator scheduling kind setter.
 * @param kind the new scheduling kind to use
 */
[[maybe_unused]] void autopas::OpenMPConfigurator::setKind(OpenMPKindOption kind) { _kind = kind; }

/**
 * Tells whether the scheduling chunk size should be overwritten.
 * @return whether the scheduling chunk size should be overwritten
 */
[[maybe_unused]] [[nodiscard]] bool autopas::OpenMPConfigurator::overrideChunkSize() const { return _kind >= 1; }

/**
 * Tells whether the scheduling kind is a manual LB4OMP scheduling technique.
 * @return whether the scheduling kind is a manual LB4OMP scheduling technique
 */
// NOLINTNEXTLINE: the function can only be static if not AUTOPAS_USE_LB4OMP.
[[maybe_unused]] [[nodiscard]] bool autopas::OpenMPConfigurator::manualSchedulingTechnique() const {
#ifdef AUTOPAS_USE_LB4OMP
  switch (_kind) {
    case OpenMPKindOption::omp_auto:
    case OpenMPKindOption::omp_dynamic:
    case OpenMPKindOption::omp_guided:
    case OpenMPKindOption::omp_runtime:
    case OpenMPKindOption::omp_static:
    case OpenMPKindOption::auto4omp_randomsel:
    case OpenMPKindOption::auto4omp_exhaustivesel:
    case OpenMPKindOption::auto4omp_expertsel:
      return false;
    default:
      return true;
  }
#else
  return false;
#endif
};
}  // namespace autopas