/**
 * @file PsVLC08CellHandler.h
 * @author Lars Doll
 * @date 20.12.2025
 */

#pragma once

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandlerUtility.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The base step processBaseCell() computes one set of pairwise interactions
 * between two cells for each spatial direction based on the baseIndex.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class PsVLC08CellHandler {
 public:
  /**
   * Constructor of the PsVLC08CellHandler.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit PsVLC08CellHandler(PairwiseFunctor *pairwiseFunctor, const std::array<unsigned long, 3> &cellsPerDimension,
                              double interactionLength, const std::array<double, 3> &cellLength,
                              const std::array<unsigned long, 3> &overlap, DataLayoutOption dataLayout, bool useNewton3)
      : _cellPairOffsets{LCC08CellHandlerUtility::computePairwiseCellOffsetsC08<
            LCC08CellHandlerUtility::C08OffsetMode::c08CellPairsSorting>(cellsPerDimension, cellLength,
                                                                         interactionLength)},
        _overlap(overlap),
        _dataLayout(dataLayout),
        _useNewton3(useNewton3),
        _cellFunctor(pairwiseFunctor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/,
                     dataLayout, useNewton3),
        _interactionLength(interactionLength),
        _cellLength(cellLength) {}

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param baseIndex Index respective to which box is constructed.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long baseIndex);

  /**
   * Sets the orientationList.
   * @param list
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell>>> &list);

 protected:
  /**
   * Pair sets for processBaseCell().
   * Values are: offset of first cell, offset of second cell, sorting direction.
   */
  std::vector<LCC08CellHandlerUtility::OffsetPairSorting> _cellPairOffsets;

  /**
   * Overlap of interacting cells. Array allows asymmetric cell sizes.
   */
  const std::array<unsigned long, 3> _overlap;

  /**
   * The datalayout to be used.
   */
  DataLayoutOption _dataLayout;

  /**
   * If newton3 should be used or not.
   */
  bool _useNewton3;

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<ParticleCell, PairwiseFunctor,
                            /*bidirectional*/ true>
      _cellFunctor;

  /**
   * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;
};

template <class ParticleCell, class PairwiseFunctor>
void PsVLC08CellHandler<ParticleCell, PairwiseFunctor>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell>>> &list) {
  _cellFunctor.setOrientationList(list);
}

template <class ParticleCell, class PairwiseFunctor>
inline void PsVLC08CellHandler<ParticleCell, PairwiseFunctor>::processBaseCell(std::vector<ParticleCell> &cells,
                                                                               unsigned long baseIndex) {
  for (auto const &[offset1, offset2, r] : _cellPairOffsets) {
    const unsigned long cellIndex1 = baseIndex + offset1;
    const unsigned long cellIndex2 = baseIndex + offset2;

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cellIndex1);
    } else {
      this->_cellFunctor.processCellPair(cellIndex1, cellIndex2, r);
    }
  }
}

}  // namespace autopas
