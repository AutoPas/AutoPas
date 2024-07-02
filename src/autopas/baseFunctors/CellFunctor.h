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

namespace autopas::internal {
/**
 * A cell functor. This functor is built from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles.
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell, class ParticleFunctor, bool bidirectional = true>
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
  explicit CellFunctor(ParticleFunctor *f, const double sortingCutoff, DataLayoutOption dataLayout, bool useNewton3)
      : _functor(f), _sortingCutoff(sortingCutoff), _dataLayout(dataLayout), _useNewton3(useNewton3) {}

  /**
   * Process the interactions inside one cell.
   * @param cell All pairwise interactions of particles inside this cell are calculated.
   */
  void processCell(ParticleCell &cell);

  /**
   * Process the interactions between the particles of cell1 with particles of cell2.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2. If no parameter or {0, 0, 0} is
   * given, sorting will be disabled.
   */
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2,
                       const std::array<double, 3> &sortingDirection = {0., 0., 0.});

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

  /**
   * Set the sorting-threshold
   * If the sum of the number of particles in two cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in two cells from which sorting should be enabled.
   */
  void setSortingThreshold(size_t sortingThreshold);

 private:
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
  void processCellAoS(ParticleCell &cell);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _sortingCutoff;

  /**
   * Min. number of particles to start sorting. This is the sum of the number of particles in two cells.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t _sortingThreshold{8};

  DataLayoutOption::Value _dataLayout;

  bool _useNewton3;
};

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::setSortingThreshold(size_t sortingThreshold) {
  _sortingThreshold = sortingThreshold;
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCell(ParticleCell &cell) {
  if ((_dataLayout == DataLayoutOption::soa and cell._particleSoABuffer.size() == 0) or
      (_dataLayout == DataLayoutOption::aos and cell.isEmpty())) {
    return;
  }

  // avoid force calculations if the cell contains only halo particles or if the cell is empty (=dummy)
  const bool cellHasOwnedParticles = toInt64(cell.getPossibleParticleOwnerships() & OwnershipState::owned);
  if (not cellHasOwnedParticles) {
    return;
  }

  switch (_dataLayout) {
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
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellPair(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
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

  switch (_dataLayout) {
    case DataLayoutOption::aos:
      if (_useNewton3) {
        processCellPairAoSN3(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoSNoN3(cell1, cell2, sortingDirection);
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
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellAoS(ParticleCell &cell) {
  // helper function
  const auto interactParticles = [&](auto &p1, auto &p2) {
    if (_useNewton3) {
      _functor->AoSFunctor(p1, p2, true);
    } else {
      _functor->AoSFunctor(p1, p2, false);
      _functor->AoSFunctor(p2, p1, false);
    }
  };

  if (cell.size() > _sortingThreshold) {
    SortedCellView<ParticleCell> cellSorted(cell, utils::ArrayMath::normalize(cell.getCellLength()));

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
  } else {
    for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
      auto p2Ptr = p1Ptr;
      ++p2Ptr;
      for (; p2Ptr != cell.end(); ++p2Ptr) {
        interactParticles(*p1Ptr, *p2Ptr);
      }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if ((cell1.size() + cell2.size() > _sortingThreshold) and (sortingDirection != std::array<double, 3>{0., 0., 0.})) {
    SortedCellView<ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell> cell2Sorted(cell2, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        _functor->AoSFunctor(*p1Ptr, *p2Ptr, true);
      }
    }
  } else {
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        _functor->AoSFunctor(p1, p2, true);
      }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellPairAoSNoN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  // helper function
  const auto interactParticlesNoN3 = [&](auto &p1, auto &p2) {
    _functor->AoSFunctor(p1, p2, false);
    if constexpr (bidirectional) {
      _functor->AoSFunctor(p2, p1, false);
    }
  };

  if ((cell1.size() + cell2.size() > _sortingThreshold) and (sortingDirection != std::array<double, 3>{0., 0., 0.})) {
    SortedCellView<ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<ParticleCell> cell2Sorted(cell2, sortingDirection);

    for (auto &[p1Projection, p1Ptr] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2Ptr] : cell2Sorted._particles) {
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        interactParticlesNoN3(*p1Ptr, *p2Ptr);
      }
    }
  } else {
    for (auto &p1 : cell1) {
      for (auto &p2 : cell2) {
        interactParticlesNoN3(p1, p2);
      }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoAN3(ParticleCell &cell1,
                                                                                     ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellPairSoANoN3(ParticleCell &cell1,
                                                                                       ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if constexpr (bidirectional) {
    _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
  }
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellSoAN3(ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class ParticleCell, class ParticleFunctor, bool bidirectional>
void CellFunctor<ParticleCell, ParticleFunctor, bidirectional>::processCellSoANoN3(ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}
}  // namespace autopas::internal
