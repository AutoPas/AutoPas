/**
 * @file OpenMPConfigurator.cpp
 * @author MehdiHachicha
 * @date 13.03.2024
 */

#include "autopas/utils/OpenMPConfigurator.h"

namespace autopas {
/**
 * OpenMP default chunk size.
 * md-flexible: set via command-line option --openmp-chunk-size <int>
 */
int openMPDefaultChunkSize = 0;

/**
 * Global OpenMP chunk size for AutoPas auto-tuning.
 */
[[maybe_unused]] int openMPTunedChunkSize = 1;

/**
 * OpenMP default chunk size.
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
    case OpenMPKindOption::auto4omp_binarySearch:
      return 4;
    case OpenMPKindOption::auto4omp_expertsel:
      return 5;
    default:
      return _chunkSize;
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
#ifdef AUTOPAS_USE_LB4OMP // LB4OMP's scheduling techniques.
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
    case OpenMPKindOption::lb4omp_viss:
      return omp_sched_viss;
    case OpenMPKindOption::lb4omp_rnd:
      return omp_sched_rnd;
    case OpenMPKindOption::lb4omp_trapezoidal:
      return omp_sched_trapezoidal;
    case OpenMPKindOption::lb4omp_fac2a:
      return omp_sched_fac2a;
    case OpenMPKindOption::lb4omp_static_steal:
      return omp_sched_static_steal;
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
    default:
      return omp_sched_auto; // Standard auto and Auto4OMP's selection methods.
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
[[maybe_unused]] [[nodiscard]] bool autopas::OpenMPConfigurator::standard() const {
  switch (_kind) {
    case OpenMPKindOption::omp_auto:
    case OpenMPKindOption::omp_dynamic:
    case OpenMPKindOption::omp_guided:
    case OpenMPKindOption::omp_runtime:
    case OpenMPKindOption::omp_static:
      return true;
    default:
      return false;
  }
};
}  // namespace autopas