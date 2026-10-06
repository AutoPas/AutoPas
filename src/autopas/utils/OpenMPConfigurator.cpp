/**
 * @file OpenMPConfigurator.cpp
 * @author MehdiHachicha
 * @date 13.03.2024
 */

#include "autopas/utils/OpenMPConfigurator.h"

#include "autopas/utils/ExceptionHandler.h"

namespace autopas {
/**
 * OpenMP default chunk size.
 * md-flexible: set via command-line option --openmp-chunk-size int
 */
int openMPDefaultChunkSize = 0;

/**
 * OpenMP default chunk size.
 * md-flexible: set via command-line option --openmp-kind kind
 */
OpenMPKindOption openMPDefaultKind = OpenMPKindOption::omp_static;

/**
 * OpenMP configurator default constructor.
 */
[[maybe_unused]] OpenMPConfigurator::OpenMPConfigurator() = default;

/**
 * OpenMP configurator constructor.
 * @param kind the OpenMP scheduling kind, LB4OMP scheduling technique, or Auto4OMP selection method to use
 * @param chunkSize the OpenMP scheduling chunk size to use
 */
[[maybe_unused]] OpenMPConfigurator::OpenMPConfigurator(OpenMPKindOption kind, int chunkSize) {
  setKind(kind);
  setChunkSize(chunkSize);
}

/**
 * AutoPas OpenMP configurator chunk size getter.
 * @return the current OpenMP chunk size
 */
[[maybe_unused]] [[nodiscard]] int OpenMPConfigurator::getChunkSize() const { return _chunkSize; }

/**
 * OpenMP chunk size getter for setting OpenMP's scheduling runtime variables.
 * @return the current OpenMP chunk size, directly usable in OpenMP's schedule setter
 */
[[maybe_unused]] [[nodiscard]] int OpenMPConfigurator::getOMPChunkSize() const {
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
[[maybe_unused]] void OpenMPConfigurator::setChunkSize(int chunkSize) { _chunkSize = chunkSize; }

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
    case OpenMPKindOption::lb4omp_trapezoidal:
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
    default:
      return omp_sched_auto;  // Standard auto and Auto4OMP's selection methods.
  }
}

/**
 * AutoPas OpenMP configurator scheduling kind setter.
 * @param kind the new scheduling kind to use
 */
[[maybe_unused]] void OpenMPConfigurator::setKind(OpenMPKindOption kind) { _kind = kind; }

/**
 * Tells whether the scheduling chunk size should be overwritten.
 * @return whether the scheduling chunk size should be overwritten
 */
[[maybe_unused]] [[nodiscard]] bool OpenMPConfigurator::overrideChunkSize() const { return _kind >= 1; }
}  // namespace autopas