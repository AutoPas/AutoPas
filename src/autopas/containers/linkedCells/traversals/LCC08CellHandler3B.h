/**
* @file LCC08CellHandler3B.h
* @author N. Deng
* @date 20.10.2023
 */

#pragma once

#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/utils/ThreeDimensionalMapping.h"


namespace autopas {

/**
* This class provides the base for traversals using the c08 base step.
*
* The base step processBaseCell() computes one set of triwise interactions
* between three cells for each spatial direction based on the baseIndex.
* After executing the base step on all cells all triwise interactions for
* all cells are done.
*
* @tparam ParticleCell the type of cells
* @tparam PairwiseFunctor The functor that defines the interaction of three particles.
* @tparam useSoA
* @tparam useNewton3
 */
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC08CellHandler3B {
 public:
  /**
  * Constructor of the LCC08CellHandler3B.
  * @param pairwiseFunctor The functor that defines the interaction of two particles.
  * @param cellsPerDimension The number of cells per dimension.
  * @param interactionLength Interaction length (cutoff + skin).
  * @param cellLength cell length.
  * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
  * in that case the interactionLength is needed!
   */
  explicit LCC08CellHandler3B(Functor *functor, const std::array<unsigned long, 3> &cellsPerDimension,
                              const double interactionLength, const std::array<double, 3> &cellLength,
                              const std::array<unsigned long, 3> &overlap)
      : _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/),
        _cellOffsets{},
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap) {
    computeOffsets(cellsPerDimension);
  }

  /**
   * Computes all interactions
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long baseIndex);

 protected:

  /**
    * Combination of triplets for processBaseCell().
   */
  std::vector<std::tuple<long, long, std::array<double, 3>>> _cellOffsets;

  /**
   * Computes triplets used in processBaseCell()
   * @param cellsPerDimension
   */
  void computeOffsets(const std::array<unsigned long, 3> &cellsPerDimension);

  /**
  * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<unsigned long, 3> _overlap;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between three cells.
   */
  internal::CellFunctor3B<typename ParticleCell::ParticleType, ParticleCell, Functor, dataLayout, false, false>
      _cellFunctor;

  Functor *_functor;

  /**
  * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
  * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08CellHandler3B<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long baseIndex) {
  for (auto const &[offset1, offset2, r] : _cellOffsets) {
    const unsigned long cellIndex1 = baseIndex + offset1;
    const unsigned long cellIndex2 = baseIndex + offset2;

    ParticleCell &cell1 = cells[cellIndex1];
    ParticleCell &cell2 = cells[cellIndex2];

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
  }
}

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3>::computeOffsets(
    const std::array<unsigned long, 3> &cellsPerDimension) {

}

}  // namespace autopas
