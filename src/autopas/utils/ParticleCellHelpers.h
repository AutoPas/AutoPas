/**
 * @file ParticleCellHelpers.h
 * @author seckler
 * @date 17.06.19
 */

#pragma once

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
  static bool checkParticleInCellAndUpdate(CellType & cellI, ParticleType & particleI) {
    for (auto iterator = cellI.begin(); iterator.isValid(); ++iterator) {
      if (iterator->getID() == particleI.getID()) {
        *iterator = particleI;
        return true;
      }
    }
    return false;
  }
}  // namespace internalinline

}  // namespace autopas
