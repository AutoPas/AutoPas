/**
 * @file markParticleAsDeleted.h
 * @author seckler
 * @date 27.07.20
 */

#pragma once

namespace autopas::internal {

/**
 * Marks a particle as deleted.
 *
 * This function is able to call the private markAsDeleted() function of the particle via a friend declaration.
 * The reason for this is that it should be as hard as possible to call Particle::markAsDeleted() from outside of
 * AutoPas, hence this free function is in the namespace internal.
 *
 * @tparam Particle
 * @param p
 * @note: This function should not be used from outside of AutoPas. Instead, use AutoPas::deleteParticle(iterator).
 */
template <typename Particle>
void markParticleAsDeleted(Particle &p) {
  p.markAsDeleted();
}

}  // namespace autopas::internal
