/**
 * @file OpenMPConfigurator.h
 * @author MehdiHachicha
 * @date 12.03.2024
 */

#pragma once

#include <cstddef>
#include <set>

#include "autopas/options/OpenMPKindOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * OpenMP default chunk size.
 * md-flexible: set via command-line option --openmp-chunk-size <int>
 */
extern int openMPDefaultChunkSize;

/**
 * Global OpenMP chunk size for AutoPas auto-tuning.
 */
[[maybe_unused]] extern int openMPTunedChunkSize;


/**
 * OpenMP default scheduling kind.
 * md-flexible: set via command-line option --openmp-kind <OpenMPKindOption>
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
   * Tells whether the scheduling kind is a standard OpenMP kind.
   * @return whether the scheduling kind is a standard OpenMP kind
   */
  [[maybe_unused]] [[nodiscard]] bool standard() const;
};  // class OpenMPConfigurator

/**
 * Sets OpenMP's runtime schedule from a given OpenMP configurator.
 * schedule(runtime) will then use them for the traversal in the concerned calling thread.
 * @param ompConfig the OpenMP configurator
 */
inline void autopas_set_schedule(autopas::OpenMPConfigurator ompConfig) {
  // If the configurator is set to omp_runtime, users are assumed to have set OMP_SCHEDULE manually.
  if (ompConfig.getKind() == OpenMPKindOption::omp_runtime) return;

  if (ompConfig.standard()) {
    autopas_set_schedule(ompConfig.getOMPKind(), ompConfig.getOMPChunkSize());
  } else {
    autopas_auto4omp_set_schedule(ompConfig.getOMPKind(), ompConfig.getOMPChunkSize());
  }
}  // void autopas_set_schedule
}  // namespace autopas

/*
 * Sources:
 * [1] https://www.computer.org/csdl/journal/td/2022/04/09524500/1wpqIcNI6YM
 * [2] https://ieeexplore.ieee.org/document/9825675
 */
