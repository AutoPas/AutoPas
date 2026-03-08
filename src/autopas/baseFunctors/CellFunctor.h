/**
 * @file CellFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
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
 * two cells of particles. Supports PsVL.
 * @tparam ParticleCell_T
 * @tparam ParticleFunctor_T the functor which
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional = true>
class CellFunctor {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This parameter indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   * @param dataLayout The data layout to be used.
   * @param useNewton3 Parameter to specify whether newton3 is used or not.
   */
  explicit CellFunctor(ParticleFunctor_T *f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All pairwise interactions of particles inside this cell are calculated.
   */
  void processCell(ParticleCell_T &cell);

  /**
   * Process the interactions inside one cell for PseudoVerletLists.
   * @param cellIndex All pairwise interactions of particles inside this cell are calculated.
   */
  void processCellPsVL(unsigned long cellIndex);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell_T &cell1, ParticleCell_T &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});
  /**
   * Process the interactions between the particles of cell1 with particles of cell2 for PseudoVerletLists.
   * @param cell1Index
   * @param cell2Index
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairPsVL(unsigned long cell1Index, unsigned long cell2Index,
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
  [[nodiscard]] static bool getBidirectional() { return bidirectional; }

  /**
   * Set the sorting-threshold
   * If the sum of the number of particles in two cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in two cells from which sorting should be enabled.
   */
  void setSortingThreshold(size_t sortingThreshold);

  /**
   * Sets the orientationList used by PseudoVerletLists.
   * @param list the orientationList
   */
  void setOrientationList(std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list);

 private:
  /**
   * Utilizes the SortedCellView to reduce the calculations needed by processCellAoS in certain cases.
   * Used for PseudoVerletLists or if the number of Particles in a cell is above the sortingThreshold.
   * @param cellSorted SortedCellView of the cell
   */
  void processCellAoSSorted(SortedCellView<ParticleCell_T> &cellSorted);

  /**
   * Applies the functor to all particle pairs exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of _useNewton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or
   * not:
   * - if _useNewton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if _useNewton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), passing
   * newton3=false.
   * @param cell
   */
  void processCellAoS(ParticleCell_T &cell);

  /**
   * Applies the functor to all particle pairs exploiting newtons third law of motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of _useNewton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or
   * not:
   * - if _useNewton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if _useNewton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), if both are non halo
   * particles, passing newton3=false. For PseudoVerletLists.
   * @param cellIndex
   */
  void processCellAoSPsVL(unsigned long cellIndex);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   * @tparam newton3 determines if Newton3 is used or not. The version of this function actually run must match
   * _useNewton3.
   */
  template <bool newton3>
  void processCellPairAoS(ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2.
   * For PseudoVerletLists.
   * @param cell1Index
   * @param cell2Index
   * @param directionIndex Normalized vector connecting centers of cell1 and cell2.
   * @tparam newton3 determines if Newton3 is used or not. The version of this function actually run must match
   * _useNewton3.
   */
  template <bool newton3>
  void processCellPairAoSPsVL(unsigned long cell1Index, unsigned long cell2Index, signed long directionIndex);

  /**
   * Utilizes the SortedCellViews to reduce the calculations needed by processCellPairAoS.
   * Used for PseudoVerletLists or for linked cells if the number of particles in the cells are above the
   * sortingThreshold.
   * @tparam newton3 determines if Newton3 is used or not. The version of this function actually run must match
   * _useNewton3.
   * @param cell1Sorted SortedCellView of cell1
   * @param cell2Sorted SortedCellView of cell2
   */
  template <bool newton3>
  void processCellPairAoSSorted(SortedCellView<ParticleCell_T> &cell1Sorted,
                                SortedCellView<ParticleCell_T> &cell2Sorted);

  /**
   * SoA version of processCellPair if newton3 is used.
   * @param cell1
   * @param cell2
   */
  void processCellPairSoAN3(ParticleCell_T &cell1, ParticleCell_T &cell2);

  /**
   * SoA version of processCellPair if newton3 is not used.
   * @param cell1
   * @param cell2
   */
  void processCellPairSoANoN3(ParticleCell_T &cell1, ParticleCell_T &cell2);

  /**
   * SoA version of processCell if newton3 is used.
   * @param cell
   */
  void processCellSoAN3(ParticleCell_T &cell);

  /**
   * SoA version of processCell if newton3 is not used.
   * @param cell
   */
  void processCellSoANoN3(ParticleCell_T &cell);

  /**
   * Calculates the interactions of two particles.
   * @tparam newton3 determines if the newton3 optimization is used. The version of this function actually run must
   * match _useNewton3.
   * @tparam oneCell determines if the particles interact within one cell
   * @param p1 particle in cell1
   * @param p2 particle in cell2
   */
  template <bool newton3, bool oneCell>
  void interactParticles(auto &p1, auto &p2);

  /**
   * Utilizes the threeToOneD mapping to calculate the directionIndex relative to the center of a 3x3 cube of cells.
   * For PseudoVerletLists.
   * @param sortingDirection
   * @return directionIndex
   */
  static signed long getDirectionIndex(const std::array<double, 3> &sortingDirection);

  /**
   * Flips negative directionIndices for PsVL as only positive directions are stored in _orientationList.
   * For PseudoVerletLists.
   * @param directionIndex
   * @return flipped directionIndex
   */
  static signed long flipDirectionIndex(signed long directionIndex);

  /**
   * The functor used for particle interactions.
   */
  ParticleFunctor_T *_functor;

  /**
   * This parameter indicates the maximal distance the sorted particles are to interact.
   */
  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting. This is the sum of the number of particles in two cells.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};

  /**
   * This parameter determines which data layout is used.
   */
  DataLayoutOption _dataLayout;

  /**
   * This parameter determines if the newton3 optimization is used.
   */
  bool _useNewton3;

  /**
   * Orientation List: For each cell, 13 sortedCellViews are stored, each of which sorts in the direction of the
   * neighboring cell. Used by PseudoVerletLists.
   */
  std::vector<std::vector<SortedCellView<ParticleCell_T>>> *_orientationList = nullptr;
};

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
signed long CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::flipDirectionIndex(
    signed long directionIndex) {
  return -(directionIndex + 2);
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
signed long CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::getDirectionIndex(
    const std::array<double, 3> &sortingDirection) {
  const auto &[x, y, z] = sortingDirection;
  // each direction get values in {-1,0,1}, according to the sign of the direction in sortingDirection.
  int xDirection = (x > 0) - (x < 0);
  int yDirection = (y > 0) - (y < 0);
  int zDirection = (z > 0) - (z < 0);

  constexpr std::array<int, 3> dims{3, 3, 3};
  // subtract an offset of 14 so the direction with index 0 is the first vector in positive direction.
  return autopas::utils::ThreeDimensionalMapping::threeToOneD(xDirection + 1, yDirection + 1, zDirection + 1, dims) -
         14;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setOrientationList(
    std::vector<std::vector<SortedCellView<ParticleCell_T>>> &list) {
  _orientationList = &list;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::setSortingThreshold(size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCell(ParticleCell_T &cell) {
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
      processCellAoS(cell);
      break;
    case DataLayoutOption::soa:
      if (_useNewton3) {
        processCellSoAN3(cell);
      } else {
        processCellSoANoN3(cell);
      }
      break;
    default:
      utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                         _dataLayout);
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPsVL(unsigned long cellIndex) {
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      processCellAoSPsVL(cellIndex);
      break;
    case DataLayoutOption::soa:
      utils::ExceptionHandler::exception("PsVL currently don't support SoA!");
      break;
    default:
      utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                         _dataLayout);
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPair(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
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

  // (Explicit) static cast required for Apple Clang (last tested version: 15.0.0)
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      if (_useNewton3) {
        processCellPairAoS<true>(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoS<false>(cell1, cell2, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if (_useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
      break;
    default:
      utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                         _dataLayout);
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairPsVL(
    unsigned long cell1Index, unsigned long cell2Index, const std::array<double, 3> &sortingDirection) {
  signed long directionIndex = getDirectionIndex(sortingDirection);

  // (Explicit) static cast required for Apple Clang (last tested version: 15.0.0)
  switch (static_cast<DataLayoutOption::Value>(_dataLayout)) {
    case DataLayoutOption::aos:
      if (_useNewton3) {
        processCellPairAoSPsVL<true>(cell1Index, cell2Index, directionIndex);
      } else {
        processCellPairAoSPsVL<false>(cell1Index, cell2Index, directionIndex);
      }
      break;
    case DataLayoutOption::soa:
      utils::ExceptionHandler::exception("PsVL currently don't support SoA!");
      break;
    default:
      utils::ExceptionHandler::exception("CellFunctor only supports AoS or SoA datalayouts. Data layout used: {}",
                                         _dataLayout);
      break;
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoSSorted(
    SortedCellView<ParticleCell_T> &cellSorted) {
  for (auto cellIter1 = cellSorted._particles.begin(); cellIter1 != cellSorted._particles.end(); ++cellIter1) {
    auto &[p1Projection, p1Ptr] = *cellIter1;
    // start inner loop ahead of the outer loop
    for (auto cellIter2 = std::next(cellIter1); cellIter2 != cellSorted._particles.end(); ++cellIter2) {
      auto &[p2Projection, p2Ptr] = *cellIter2;
      if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
        break;
      }
      if (_useNewton3) {
        interactParticles<true, true>(*p1Ptr, *p2Ptr);
      } else {
        interactParticles<false, true>(*p1Ptr, *p2Ptr);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoS(ParticleCell_T &cell) {
  if (cell.size() > _sortingThreshold) {
    SortedCellView<ParticleCell_T> cellSorted(cell, utils::ArrayMath::normalize(cell.getCellLength()));
    processCellAoSSorted(cellSorted);
  } else {
    for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
      auto p2Ptr = p1Ptr;
      ++p2Ptr;
      for (; p2Ptr != cell.end(); ++p2Ptr) {
        if (_useNewton3) {
          interactParticles<true, true>(*p1Ptr, *p2Ptr);
        } else {
          interactParticles<false, true>(*p1Ptr, *p2Ptr);
        }
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellAoSPsVL(unsigned long cellIndex) {
  auto &cellSorted = (*_orientationList)[cellIndex][6];  // Index 6 for the cell diagonal (direction {1,1,1})
  processCellAoSSorted(cellSorted);
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoS(
    ParticleCell_T &cell1, ParticleCell_T &cell2, const std::array<double, 3> &sortingDirection) {
  if ((cell1.size() + cell2.size() > _sortingThreshold) and (sortingDirection != std::array<double, 3>{0., 0., 0.})) {
    SortedCellView<ParticleCell_T> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell_T> cell2Sorted(cell2, sortingDirection);
    processCellPairAoSSorted<newton3>(cell1Sorted, cell2Sorted);
  } else {
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        interactParticles<newton3, false>(p1, p2);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoSPsVL(unsigned long cell1Index,
                                                                                           unsigned long cell2Index,
                                                                                           signed long directionIndex) {
  bool reverseCell2View = false;
  if (directionIndex < 0) {
    reverseCell2View = true;
    directionIndex = flipDirectionIndex(directionIndex);
  }
  auto &cell1Sorted = (*_orientationList)[cell1Index][directionIndex];
  auto &cell2Sorted = (*_orientationList)[cell2Index][directionIndex];

  if (reverseCell2View == false) {
    processCellPairAoSSorted<newton3>(cell1Sorted, cell2Sorted);
  } else {
    for (auto &[p1Proj, p1Ptr] : cell1Sorted._particles) {
      for (auto it = cell2Sorted._particles.rbegin(); it != cell2Sorted._particles.rend(); ++it) {
        const auto &[p2Proj, p2Ptr] = *it;
        if (std::abs(p1Proj - p2Proj) > _sortingCutoff) {
          break;
        }
        interactParticles<newton3, false>(*p1Ptr, *p2Ptr);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3, bool oneCell>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::interactParticles(auto &p1, auto &p2) {
  if constexpr (newton3) {
    _functor->AoSFunctor(p1, p2, true);
  } else {
    if constexpr (oneCell) {
      if (not p1.isHalo()) {
        _functor->AoSFunctor(p1, p2, false);
      }
      if (not p2.isHalo()) {
        _functor->AoSFunctor(p2, p1, false);
      }
    } else {
      _functor->AoSFunctor(p1, p2, false);
      if constexpr (bidirectional) {
        _functor->AoSFunctor(p2, p1, false);
      }
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
template <bool newton3>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairAoSSorted(
    SortedCellView<ParticleCell_T> &cell1Sorted, SortedCellView<ParticleCell_T> &cell2Sorted) {
  for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
    for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
      if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
        break;
      }
      interactParticles<newton3, false>(*p1Ptr, *p2Ptr);
    }
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoAN3(ParticleCell_T &cell1,
                                                                                         ParticleCell_T &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellPairSoANoN3(ParticleCell_T &cell1,
                                                                                           ParticleCell_T &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if constexpr (bidirectional) {
    _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
  }
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellSoAN3(ParticleCell_T &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class ParticleCell_T, class ParticleFunctor_T, bool bidirectional>
void CellFunctor<ParticleCell_T, ParticleFunctor_T, bidirectional>::processCellSoANoN3(ParticleCell_T &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}
}  // namespace autopas::internal