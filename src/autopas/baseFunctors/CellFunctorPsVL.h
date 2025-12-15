/**
 * @file CellFunctorPsVL.h
 *
 * @date 06.12.2025
 * @author Lars Doll
 */

#pragma once

#include "autopas/cells/SortedCellView.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles using the PseudoVerletList algorithm.
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell, class ParticleFunctor, bool bidirectional = true>
class CellFunctorPsVL {
 public:
  /**
   * The constructor of CellFunctorPsVL.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This parameter normally should be the cutoff + skin.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit CellFunctorPsVL(ParticleFunctor *f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All pairwise interactions of particles inside this cell are calculated.
   * @param cellIndex
   */
  void processCell(ParticleCell &cell, unsigned long cellIndex);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1
   * @param cell1Index
   * @param cell2
   * @param cell2Index
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell &cell1, unsigned long cell1Index,
    ParticleCell &cell2, unsigned long cell2Index,
    const std::array<double, 3> &sortingDirection);

  /**
   * Getter
   * @return
   */
  [[nodiscard]] DataLayoutOption getDataLayout() const { return _dataLayout; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getNewton3() const { return _useNewton3; }

  /**
   * Getter
   * @return
   */
  [[nodiscard]] bool getBidirectional() const { return bidirectional; }

 private:

  void setOrientationLists(std::vector<std::vector<SortedCellView<typename ParticleCell::ParticleType>>> &lists);

  /**
   * Applies the functor to all particle pairs exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of _useNewton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or
   * not:
   * - if _useNewton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if _useNewton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), passing
   * newton3=false.
   * @param cell
   * @param cellIndex
   */
  void processCellAoS(ParticleCell &cell, unsigned long cellIndex);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell1Index
   * @param cell2
   * @param cell2Index
   * @param directionIndex Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSN3(ParticleCell &cell1, unsigned long cell1Index, ParticleCell &cell2, unsigned long cell2Index, signed long directionIndex);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell1Index
   * @param cell2
   * @param cell2Index
   * @param directionIndex Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSNoN3(ParticleCell &cell1, unsigned long cell1Index, ParticleCell &cell2, unsigned long cell2Index, signed long directionIndex);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  DataLayoutOption _dataLayout;

  bool _useNewton3;

  std::vector<std::vector<SortedCellView<typename ParticleCell::ParticleType>>>* _orientationList = nullptr;
};

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::setOrientationLists(
    std::vector<std::vector<SortedCellView<typename ParticleCell::ParticleType>>> &lists) {
  _orientationList = &lists;
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCell(ParticleCell &cell, unsigned long cellIndex) {
  if ((_dataLayout == DataLayoutOption::soa and cell._particleSoABuffer.size() == 0) or
    (_dataLayout == DataLayoutOption::aos and cell.isEmpty())) {
    return;
    }

  // avoid force calculations if the cell contains only halo particles or if the cell is empty (=dummy)
  const bool cellHasOwnedParticles = toInt64(cell.getPossibleParticleOwnerships() & OwnershipState::owned);
  if (not cellHasOwnedParticles) {
    return;
  }

  // (Explicit) static cast required for Apple Clang (last tested version: 15.0.0)
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      processCellAoS(cell, cellIndex);
      break;
    case DataLayoutOption::soa:
      if (_useNewton3) {
        processCellSoAN3(cell, cellIndex);
      } else {
        processCellSoANoN3(cell, cellIndex);
      }
      break;
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPair(
    ParticleCell &cell1, unsigned long cell1Index, ParticleCell &cell2, unsigned long cell2Index, const std::array<double, 3> &sortingDirection) {
  if ((_dataLayout == DataLayoutOption::soa and
       (cell1._particleSoABuffer.size() == 0 or cell2._particleSoABuffer.size() == 0)) or
      (_dataLayout == DataLayoutOption::aos and (cell1.isEmpty() or cell2.isEmpty()))) {
    return;
      }

  // avoid force calculations if both cells can not contain owned particles or if newton3==false and cell1 does not
  // contain owned particles
  const bool cell1HasOwnedParticles = toInt64(cell1.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell2HasOwnedParticles = toInt64(cell2.getPossibleParticleOwnerships() & OwnershipState::owned);

  if (((not cell1HasOwnedParticles) and (not _useNewton3) and (not bidirectional)) or
      ((not cell1HasOwnedParticles) and (not cell2HasOwnedParticles))) {
    return;
      }

  const auto& [x, y, z] = sortingDirection;
  int xInt = (x > 0) - (x < 0);
  int yInt = (y > 0) - (y < 0);
  int zInt = (z > 0) - (z < 0);

  constexpr std::array<int, 3> dims{3, 3, 3};
  signed long directionIndex = autopas::utils::ThreeDimensionalMapping::threeToOneD(xInt+1,yInt+1,zInt+1,dims) - 14;

  // (Explicit) static cast required for Apple Clang (last tested version: 15.0.0)
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      if (_useNewton3) {
        processCellPairAoSN3(cell1, cell1Index, cell2, cell2Index, directionIndex);
      } else {
        processCellPairAoSNoN3(cell1, cell1Index, cell2, cell2Index, directionIndex);
      }
      break;
    case DataLayoutOption::soa:
      if (_useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
      break;
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellAoS(ParticleCell &cell, unsigned long cellIndex) {
  // helper function
  const auto interactParticles = [&](auto &p1, auto &p2) {
    if (_useNewton3) {
      _functor->AoSFunctor(p1, p2, true);
    } else {
      if (not p1.isHalo()) {
        _functor->AoSFunctor(p1, p2, false);
      }
      if (not p2.isHalo()) {
        _functor->AoSFunctor(p2, p1, false);
      }
    }
  };

  SortedCellView<ParticleCell> cellSorted= _orientationList[cellIndex][6]; //Index 6 for the cell diagonal (direction {1,1,1})

  for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
    auto &[p1Projection, p1Ptr] = *cellIter1;
    // start inner loop ahead of the outer loop
    for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
      auto &[p2Projection, p2Ptr] = *cellIter2;
      if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
        break;
      }
      interactParticles(*p1Ptr, *p2Ptr);
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, unsigned long cell1Index, ParticleCell &cell2, unsigned long cell2Index, signed long directionIndex) {
  bool revSecondView = false;
  if (directionIndex<0) {
    revSecondView = true;
    directionIndex = -(directionIndex+2);
  }
  SortedCellView<ParticleCell> cell1Sorted = _orientationList[cell1Index][directionIndex];
  SortedCellView<ParticleCell> cell2Sorted = _orientationList[cell2Index][directionIndex];

  for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {

    if (revSecondView == false) {
      // Forward
      for (auto it = cell2Sorted._particles.begin(); it != cell2Sorted._particles.end(); ++it) {
        const auto &[p2Projection, p2Ptr] = *it;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        _functor->AoSFunctor(*p1Ptr, *p2Ptr, true);
           }
    } else {
      // Backward
      for (auto it = cell2Sorted._particles.rbegin(); it != cell2Sorted._particles.rend(); ++it) {
        const auto &[p2Projection, p2Ptr] = *it;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        _functor->AoSFunctor(*p1Ptr, *p2Ptr, true);
           }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSNoN3(
    ParticleCell &cell1, unsigned long cell1Index, ParticleCell &cell2, unsigned long cell2Index, signed long directionIndex) {
  const auto interactParticlesNoN3 = [&](auto &p1, auto &p2) {
    _functor->AoSFunctor(p1, p2, false);
    if constexpr (bidirectional) {
      _functor->AoSFunctor(p2, p1, false);
    }
  };

  bool revSecondView = false;
  if (directionIndex<0) {
    revSecondView = true;
    directionIndex = -(directionIndex+2);
  }

  SortedCellView<ParticleCell> cell1Sorted = _orientationList[cell1Index][directionIndex];
  SortedCellView<ParticleCell> cell2Sorted = _orientationList[cell2Index][directionIndex];

  for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {

    if (revSecondView == false) {
      // Forward
      for (auto it = cell2Sorted._particles.begin(); it != cell2Sorted._particles.end(); ++it) {
        const auto &[p2Projection, p2Ptr] = *it;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        interactParticlesNoN3(*p1Ptr, *p2Ptr);
           }
    } else {
      // Backward
      for (auto it = cell2Sorted._particles.rbegin(); it != cell2Sorted._particles.rend(); ++it) {
        const auto &[p2Projection, p2Ptr] = *it;
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        interactParticlesNoN3(*p1Ptr, *p2Ptr);
           }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoAN3(ParticleCell &cell1,
                                                                                     ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoANoN3(ParticleCell &cell1,
                                                                                       ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if constexpr (bidirectional) {
    _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellSoAN3(ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctorPsVL<ParticleCell, ParticleFunctor, bidirectional>::processCellSoANoN3(ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}
}  // namespace autopas::internal
