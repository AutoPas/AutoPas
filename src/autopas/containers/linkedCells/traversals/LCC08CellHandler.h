/**
 * @file LCC08CellHandler.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandlerUtility.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/checkFunctorType.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The base step processBaseCell() computes one set of pairwise interactions
 * between two cells for each spatial direction based on the baseIndex.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam Functor_T The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class Functor_T>
class LCC08CellHandler {
 public:
  /**
   * Constructor of the LCC08CellHandler.
   * @param functor The functor that defines the interaction of two particles.
   * @param cellsPerDimension The number of cells per dimension.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param overlap number of overlapping cells in each direction as result from cutoff and cellLength.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC08CellHandler(Functor_T *functor, const std::array<unsigned long, 3> &cellsPerDimension,
                            double interactionLength, const std::array<double, 3> &cellLength,
                            const std::array<unsigned long, 3> &overlap, DataLayoutOption dataLayout, bool useNewton3)
      : _cellFunctor(functor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/,
                     dataLayout, useNewton3),
        _interactionLength(interactionLength),
        _cellLength(cellLength),
        _overlap(overlap),
        _dataLayout(dataLayout),
        _useNewton3(useNewton3) {
    if constexpr (utils::isPairwiseFunctor<Functor_T>()) {
      _cellOffsets =
          LCC08CellHandlerUtility::computePairwiseCellOffsetsC08<LCC08CellHandlerUtility::C08OffsetMode::sorting>(
              cellsPerDimension, cellLength, interactionLength);
    } else if constexpr (utils::isTriwiseFunctor<Functor_T>()) {
      _cellOffsets =
          LCC08CellHandlerUtility::computeTriwiseCellOffsetsC08<LCC08CellHandlerUtility::C08OffsetMode::sorting>(
              cellsPerDimension, cellLength, interactionLength);
    } else {
      utils::ExceptionHandler::exception("LCC08CellHandler::LCC08CellHandler(): Functor is not valid.");
    }
  }

  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param baseIndex Index respective to which box is constructed.
   */
  void processBaseCell(std::vector<ParticleCell_T> &cells, unsigned long baseIndex);

  /**
   * Pairwise implementation of processBaseCell().
   * @copydoc processBaseCell()
   */
  inline void processBaseCellPairwise(std::vector<ParticleCell_T> &cells, unsigned long baseIndex);

  /**
   * Triwise implementation of processBaseCell().
   * @copydoc processBaseCell()
   */
  inline void processBaseCellTriwise(std::vector<ParticleCell_T> &cells, unsigned long baseIndex);

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) { _cellFunctor.setSortingThreshold(sortingThreshold); }

 protected:
  // CellOffset needs to store interaction pairs or triplets depending on the Functor type.
  using CellOffsetType =
      std::conditional_t<decltype(utils::isPairwiseFunctor<Functor_T>())::value,
                         LCC08CellHandlerUtility::OffsetPairSorting, LCC08CellHandlerUtility::OffsetTripletSorting>;

  /**
   * Pair sets for processBaseCell().
   * Values are: offset of first cell, offset of second cell, sorting direction.
   */
  std::vector<CellOffsetType> _cellOffsets;

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
  // CellFunctor type for either Pairwise or Triwise Functors.e
  using CellFunctorType =
      std::conditional_t<decltype(utils::isPairwiseFunctor<Functor_T>())::value,
                         internal::CellFunctor<ParticleCell_T, Functor_T, /*bidirectional*/ true>,
                         internal::CellFunctor3B<ParticleCell_T, Functor_T, /*bidirectional*/ true>>;

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctorType _cellFunctor;

  /**
   * Interaction length (cutoff + skin).
   */
  const double _interactionLength;

  /**
   * Cell length in CellBlock3D.
   */
  const std::array<double, 3> _cellLength;
};

template <class ParticleCell_T, class Functor_T>
inline void LCC08CellHandler<ParticleCell_T, Functor_T>::processBaseCell(std::vector<ParticleCell_T> &cells,
                                                                         unsigned long baseIndex) {
  if constexpr (utils::isPairwiseFunctor<Functor_T>()) {
    processBaseCellPairwise(cells, baseIndex);
  } else if constexpr (utils::isTriwiseFunctor<Functor_T>()) {
    processBaseCellTriwise(cells, baseIndex);
  } else {
    utils::ExceptionHandler::exception(
        "LCC01Traversal::processBaseCell(): Given Functor type is not of type PairwiseFunctor or TriwiseFunctor.");
  }
}

template <class ParticleCell_T, class Functor_T>
inline void LCC08CellHandler<ParticleCell_T, Functor_T>::processBaseCellPairwise(std::vector<ParticleCell_T> &cells,
                                                                                 unsigned long baseIndex) {
  for (auto const &[offset1, offset2, r] : _cellOffsets) {
    const unsigned long cellIndex1 = baseIndex + offset1;
    const unsigned long cellIndex2 = baseIndex + offset2;

    ParticleCell_T &cell1 = cells[cellIndex1];
    ParticleCell_T &cell2 = cells[cellIndex2];

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2, r);
    }
  }
}

template <class ParticleCell_T, class Functor_T>
inline void LCC08CellHandler<ParticleCell_T, Functor_T>::processBaseCellTriwise(std::vector<ParticleCell_T> &cells,
                                                                                unsigned long baseIndex) {
  for (auto const &[offset1, offset2, offset3, r] : _cellOffsets) {
    const unsigned long index1 = baseIndex + offset1;
    const unsigned long index2 = baseIndex + offset2;
    const unsigned long index3 = baseIndex + offset3;

    ParticleCell_T &cell1 = cells[index1];
    ParticleCell_T &cell2 = cells[index2];
    ParticleCell_T &cell3 = cells[index3];

    if (index1 == index2 && index1 == index3 && index2 == index3) {
      this->_cellFunctor.processCell(cell1);
    } else if (index1 == index2 && index1 != index3) {
      this->_cellFunctor.processCellPair(cell1, cell3);
    } else if (index1 != index2 && index1 == index3) {
      this->_cellFunctor.processCellPair(cell1, cell2);
    } else if (index1 != index2 && index2 == index3) {
      this->_cellFunctor.processCellPair(cell1, cell2);
    } else {
      this->_cellFunctor.processCellTriple(cell1, cell2, cell3);
    }
  }
}

}  // namespace autopas
