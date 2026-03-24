/**
 * @file HierarchicalCellBlock3D.h
 *
 * @date 21 Mar 2026
 */

#pragma once

#include "autopas/containers/CellBlock3D.h"

namespace autopas::internal {

/**
 * CellBlock variant for HierarchicalGrid containers.
 * @tparam ParticleCell Type of the handled ParticleCells.
 */
template <class ParticleCell>
class HierarchicalCellBlock3D : public CellBlock3D<ParticleCell> {
 public:
  using typename CellBlock3D<ParticleCell>::index_t;

  /**
   * Constructor of HierarchicalCellBlock3D.
   * @param vec Vector of ParticleCells that this class manages.
   * @param bMin Lower corner of the cell block.
   * @param bMax Higher corner of the cell block.
   * @param interactionLength Max. radius of interaction between particles.
   * @param cellSizeFactor Cell size factor relative to interactionLength.
   */
  HierarchicalCellBlock3D(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin,
                          const std::array<double, 3> &bMax, double interactionLength,
                          std::array<size_t, 3> cellsPerDim)
      : _cells(&vec), _boxMin(bMin), _boxMax(bMax), _interactionLength(interactionLength) {
    rebuild(vec, bMin, bMax, interactionLength, cellsPerDim);
  }

  /**
   * Deleted copy constructor.
   */
  HierarchicalCellBlock3D(const HierarchicalCellBlock3D &) = delete;

  /**
   * Deleted assignment operator.
   * @return
   */
  HierarchicalCellBlock3D &operator=(const HierarchicalCellBlock3D &) = delete;

  /**
   * Rebuild the cell block. This resizes the cell block and sets the appropriate internal variables.
   * @param vec New vector of ParticleCells to which the internal pointer is set.
   * @param bMin New lower corner of the cell block.
   * @param bMax New higher corner of the cell block.
   * @param interactionLength Max. radius of interaction between particles.
   * @param cellsPerDim Number of cells per dimension.
   */
  void rebuild(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin, const std::array<double, 3> &bMax,
               double interactionLength, std::array<size_t, 3> &cellsPerDim);
};

template <class ParticleCell>
inline void HierarchicalCellBlock3D<ParticleCell>::rebuild(std::vector<ParticleCell> &vec,
                                                           const std::array<double, 3> &bMin,
                                                           const std::array<double, 3> &bMax, double interactionLength,
                                                           std::array<size_t, 3> &cellsPerDim) {
  using namespace autopas::utils::ArrayMath::literals;
using index_t = typename CellBlock3D<ParticleCell>::index_t;                                                 
  this->_cells = &vec;
  this->_boxMin = bMin;
  this->_boxMax = bMax;
  this->_interactionLength = interactionLength;

  this->_cellsPerInteractionLength = 1;

  // compute cell length and number of cells
  index_t numCells = 1;
  for (int d = 0; d < 3; ++d) {
    const double boxLength = this->_boxMax[d] - this->_boxMin[d];

    this->_cellLength[d] = boxLength / static_cast<double>(cellsPerDim[d]);

    this->_cellLengthReciprocal[d] = static_cast<double>(cellsPerDim[d]) / boxLength;  // compute with least rounding possible

    auto tmp= static_cast<index_t>(std::ceil(this->_interactionLength / this->_cellLength[d]));

    if(tmp> this->_cellsPerInteractionLength) {
      this->_cellsPerInteractionLength = tmp;
    }
  }

  for(int d = 0; d < 3; ++d) {
    this->_cellsPerDimensionWithHalo[d] = cellsPerDim[d] + 2 * this->_cellsPerInteractionLength;

    this->_haloBoxMin[d] = this->_boxMin[d] - this->_cellsPerInteractionLength * this->_cellLength[d];
    this->_haloBoxMax[d] = this->_boxMax[d] + this->_cellsPerInteractionLength * this->_cellLength[d];

    numCells *= this->_cellsPerDimensionWithHalo[d];
  }

  AutoPasLog(TRACE, "Box Length incl Halo : {}", autopas::utils::ArrayUtils::to_string(this->_haloBoxMax - this->_haloBoxMin));
  AutoPasLog(TRACE, "Cells/Dim  incl Halo : {}", autopas::utils::ArrayUtils::to_string(this->_cellsPerDimensionWithHalo));
  AutoPasLog(TRACE, "Cell Length          : {}", autopas::utils::ArrayUtils::to_string(this->_cellLength));
  AutoPasLog(TRACE, "Interaction Length   : {}", interactionLength);

  this->_firstOwnedCellIndex = this->_cellsPerDimensionWithHalo[0] * this->_cellsPerDimensionWithHalo[1] * this->_cellsPerInteractionLength +
                         this->_cellsPerDimensionWithHalo[0] * this->_cellsPerInteractionLength + this->_cellsPerInteractionLength;
  this->_lastOwnedCellIndex = numCells - 1 - this->_firstOwnedCellIndex;

  // initialize cells
  this->_cells->resize(numCells);

  for (auto &cell : *this->_cells) {
    cell.setCellLength(this->_cellLength);
  }

  // determine the OwnershipStates, each cell can contain. This is later used in the CellFunctor to skip calculations
  for (int i = 0; i < numCells; i++) {
    const bool canHaveHalos = cellCanContainHaloParticles(i);
    const bool canHaveOwned = cellCanContainOwnedParticles(i);
    if (canHaveHalos and canHaveOwned) {
      (*this->_cells)[i].setPossibleParticleOwnerships(OwnershipState::owned | OwnershipState::halo);
    } else if (canHaveHalos) {
      (*this->_cells)[i].setPossibleParticleOwnerships(OwnershipState::halo);
    } else if (canHaveOwned) {
      (*this->_cells)[i].setPossibleParticleOwnerships(OwnershipState::owned);
    } else {
      (*this->_cells)[i].setPossibleParticleOwnerships(OwnershipState::dummy);
    }
  }
}

}  // namespace autopas::internal
