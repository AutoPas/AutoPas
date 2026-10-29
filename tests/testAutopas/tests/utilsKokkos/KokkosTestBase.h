/**
 * @file KokkosTestBase.h
 * @author Luis Gall
 * @date 29.10.2026
 */

#pragma once

#ifdef AUTOPAS_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>

#include "AutoPasTestBase.h"

/**
 * Base fixture for all tests that need a Kokkos runtime.
 *
 * Kokkos can only be initialized once per process, and not again after Kokkos::finalize().
 * Therefore, the runtime is initialized lazily by the first Kokkos test suite and finalized once at the end of main().
 */
class KokkosTestBase : public AutoPasTestBase {
 public:
  /**
   * Initializes Kokkos if no previous test suite has done so yet.
   */
  static void SetUpTestSuite() {
    if (not Kokkos::is_initialized()) {
      Kokkos::initialize();
    }
  }
};

#endif