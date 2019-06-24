/**
 * @file ParticleCellHelpers.h
 * @author seckler
 * @date 17.06.19
 */

#pragma once
#include "autopas/utils/ArrayMath.h"

namespace autopas {

namespace internal {
/**
 * Updates a found particle within cellI to the values of particleI.
 * Checks whether a particle with the same id as particleI is within the cell
 * cellI and overwrites the particle with particleI, if it is found.
 * @param cellI
 * @param particleI
 * @tparam ParticleType
 * @tparam CellType
 * @return
 */
template <class ParticleType, class CellType>
static bool checkParticleInCellAndUpdate(CellType &cellI, ParticleType &particleI) {
  for (auto iterator = cellI.begin(); iterator.isValid(); ++iterator) {
    if (iterator->getID() == particleI.getID()) {
      *iterator = particleI;
      return true;
    }
  }
  return false;
}

/**
 * Same as checkParticleInCellAndUpdate(CellType, ParticleType), but additionally checks whether the particle is close
 * to the other particle:
 * @copydoc checkParticleInCellAndUpdate()
 * @param absError maximal distance the previous particle is allowed to be away from the new particle.
 * @note this version is useful, if there might be more than one particle with the same id in the same cell.
 */
template <class ParticleType, class CellType>
static bool checkParticleInCellAndUpdateNearPosition(CellType &cellI, ParticleType &particleI, double absError) {
  for (auto iterator = cellI.begin(); iterator.isValid(); ++iterator) {
    if (iterator->getID() == particleI.getID()) {
      auto distanceVec = autopas::ArrayMath::sub(iterator->getR(), particleI.getR());
      auto distanceSqr = autopas::ArrayMath::dot(distanceVec, distanceVec);
      if (distanceSqr < absError * absError) {
        *iterator = particleI;
        // found the particle, returning.
        return true;
      }
    }
  }
  return false;
}
}  // namespace internal
}  // namespace autopas
