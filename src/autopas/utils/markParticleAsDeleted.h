/**
 * @file markParticleAsDeleted.h
 * @author seckler
 * @date 27.07.20
 */

#pragma once

namespace autopas::internal {

/**
 * Marks a particle as deleted.
 * @tparam Particle
 * @param p
 * @note: This function should not be used from outside of AutoPas. Instead, use AutoPas::deleteParticle(iterator).
 */
template <typename Particle>
void markParticleAsDeleted(Particle &p) {
  p.markAsDeleted();
}

}  // namespace autopas::internal
