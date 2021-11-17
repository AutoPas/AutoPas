/**
 * @file markParticleAsDeleted.h
 * @author seckler
 * @date 27.07.20
 */

#pragma once
#include <Kokkos_Core.hpp>

namespace autopas::internal {

/**
 * Marks a particle as deleted.
 * @tparam Particle
 * @param p
 * @note: This function should not be used from outside of AutoPas. Instead, use AutoPas::deleteParticle(iterator).
 */
template <typename Particle>
KOKKOS_INLINE_FUNCTION void markParticleAsDeleted(Particle &p) {
  p.markAsDeleted();
}

}  // namespace autopas::internal