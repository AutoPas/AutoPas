/**
 * @file ParticleCellHelpers.h
 * @author seckler
 * @date 17.06.19
 */

#pragma once
#include "autopas/utils/ArrayMath.h"

namespace autopas::internal {
/**
 * Updates a found particle within cellI to the values of particleI.
 * Checks whether a particle with the same id as particleI is within the cell
 * cellI and overwrites the particle with particleI, if it is found.
 * @param cell
 * @param particle
 * @tparam ParticleType
 * @tparam CellType
 * @return true if the particle was updated, false otherwise.
 */
template <class ParticleType, class CellType>
static bool checkParticleInCellAndUpdateByID(CellType &cell, const ParticleType &particle) {
  for (auto iterator = cell.begin(); iterator.isValid(); ++iterator) {
    if (iterator->getID() == particle.getID()) {
      *iterator = particle;
      return true;
    }
  }
  return false;
}

/**
 * Same as checkParticleInCellAndUpdateByID(CellType, ParticleType), but additionally checks whether the particle is
 * close to the other particle:
 * @copydoc checkParticleInCellAndUpdateByID()
 * @param absError maximal distance the previous particle is allowed to be away from the new particle.
 * @note This version is useful, if there might be more than one particle with the same id in the same cell.
 */
template <class ParticleType, class CellType>
static bool checkParticleInCellAndUpdateByIDAndPosition(CellType &cell, const ParticleType &particle, double absError) {
  using namespace autopas::utils::ArrayMath::literals;
  for (auto iterator = cell.begin(); iterator.isValid(); ++iterator) {
    if (iterator->getID() == particle.getID()) {
      auto distanceVec = iterator->getR() - particle.getR();
      auto distanceSqr = autopas::utils::ArrayMath::dot(distanceVec, distanceVec);
      if (distanceSqr < absError * absError) {
        *iterator = particle;
        // found the particle, returning.
        return true;
      }
    }
  }
  return false;
}
}  // namespace autopas::internal
