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
 * cellI and close to particleI, and overwrites the particle with particleI, if it is found.
 * @param cell
 * @param particle
 * @param absError maximal distance the previous particle is allowed to be away from the new particle.
 * @tparam CellType
 * @return true if the particle was updated, false otherwise.
 * @note The position is checked, because there might be more than one particle with the same id in the same cell or in
 * neighboring cells.
 * @note The cell is not locked. If iterating the cell is not thread safe, the caller has to lock it.
 */
template <class CellType>
static bool checkParticleInCellAndUpdateByIDAndPosition(CellType &cell, const typename CellType::ParticleType &particle,
                                                        double absError) {
  using namespace autopas::utils::ArrayMath::literals;
  for (auto &p : cell) {
    if (p.getID() == particle.getID()) {
      auto distanceVec = p.getR() - particle.getR();
      auto distanceSqr = autopas::utils::ArrayMath::dot(distanceVec, distanceVec);
      if (distanceSqr < absError * absError) {
        p = particle;
        // found the particle, returning.
        return true;
      }
    }
  }
  return false;
}
}  // namespace autopas::internal
