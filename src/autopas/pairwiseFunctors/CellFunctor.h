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
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam dataLayout the dataLayout to be used
 * @tparam useNewton3
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3 = true, bool bidirectional = true>
class CellFunctor {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This paramater indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   */
  explicit CellFunctor(ParticleFunctor *f, const double sortingCutoff) : _functor(f), _sortingCutoff(sortingCutoff) {}

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
  DataLayoutOption::Value getDataLayout() const {
    return dataLayout;
  }

  /**
   * Getter
   * @return
   */
  bool getNewton3() const {
    return useNewton3;
  }

  /**
   * Getter
   */
  bool getBidirectional() const {
    return bidirectional;
  }

 private:
  /**
   * Applies the functor to all particle pairs exploiting newtons third law of
   * motion.
   * There is only one version of this function as newton3 is always allowed to be applied inside of a cell.
   * The value of newton3 defines whether or whether not to apply the aos version functor in a newton3 fashion or not:
   * - if newton3 is true: the aos functor will be applied once for each pair (only i,j), passing newton3=true.
   * - if newton3 is false: the aos functor will be applied twice for each pair (i,j and j,i), passing newton3=false.
   * @tparam newton3 defines whether or not to use newton3
   * @param cell
   */
  template <bool newton3>
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
   * Min. number of particles to start sorting.
   * @todo Currently, this is disabled because of https://github.com/AutoPas/AutoPas/issues/418
   */
  constexpr static unsigned long _startSorting = std::numeric_limits<unsigned long>::max();
};

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCell(
    ParticleCell &cell) {
  if ((dataLayout == DataLayoutOption::soa and cell._particleSoABuffer.size() == 0) or
      (dataLayout == DataLayoutOption::aos and cell.size() == 0)) {
    return;
  }

  // avoid force calculations if the cell contains only halo particles or if the cell is empty (=dummy)
  const bool cellHasOwnedParticles = toInt64(cell.getPossibleParticleOwnerships() & OwnershipState::owned);
  if (not cellHasOwnedParticles) {
    return;
  }

  switch (dataLayout) {
    case DataLayoutOption::aos:
      processCellAoS<useNewton3>(cell);
      break;
    case DataLayoutOption::soa:
      if constexpr (useNewton3) {
        processCellSoAN3(cell);
      } else {
        processCellSoANoN3(cell);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCellPair(

    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if ((dataLayout == DataLayoutOption::soa and
       (cell1._particleSoABuffer.size() == 0 and cell2._particleSoABuffer.size() == 0)) or
      (dataLayout == DataLayoutOption::aos and (cell1.size() == 0 and cell2.size() == 0))) {
    return;
  }

  // avoid force calculations if both cells can not contain owned particles or if newton3==false and cell1 does not
  // contain owned particles
  const bool cell1HasOwnedParticles = toInt64(cell1.getPossibleParticleOwnerships() & OwnershipState::owned);
  const bool cell2HasOwnedParticles = toInt64(cell2.getPossibleParticleOwnerships() & OwnershipState::owned);

  if (((not cell1HasOwnedParticles) and (not useNewton3) and (not bidirectional)) or
      ((not cell1HasOwnedParticles) and (not cell2HasOwnedParticles))) {
    return;
  }

  switch (dataLayout) {
    case DataLayoutOption::aos:
      if constexpr (useNewton3) {
        processCellPairAoSN3(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoSNoN3(cell1, cell2, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if constexpr (useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
template <bool newton3>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCellAoS(
    ParticleCell &cell) {
  // helper function
  const auto interactParticles = [&](auto &p1, auto &p2) {
    if constexpr (newton3) {
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

  if (cell.size() > _startSorting) {
    SortedCellView<Particle, ParticleCell> cellSorted(
        cell, utils::ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

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
    for (auto cellIter1 = cell.begin(); cellIter1 != cell.end(); ++cellIter1) {
      auto cellIter2 = cellIter1;
      ++cellIter2;
      for (; cellIter2 != cell.end(); ++cellIter2) {
        interactParticles(*cellIter1, *cellIter2);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &sortingDirection) {
  if (cell1.size() + cell2.size() > _startSorting and sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

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

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3,
                 bidirectional>::processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                                                        const std::array<double, 3> &sortingDirection) {
  // helper function
  const auto interactParticlesNoN3 = [&](auto &p1, auto &p2) {
    _functor->AoSFunctor(p1, p2, false);
    if constexpr (bidirectional) {
      if (p2.isOwned()) {
        _functor->AoSFunctor(p2, p1, false);
      }
    }
  };

  if (cell1.size() + cell2.size() > _startSorting and sortingDirection != std::array<double, 3>{0., 0., 0.}) {
    SortedCellView<Particle, ParticleCell> cell1Sorted(cell1, sortingDirection);
    SortedCellView<Particle, ParticleCell> cell2Sorted(cell2, sortingDirection);

    for (auto &[p1Projection, p1] : cell1Sorted._particles) {
      for (auto &[p2Projection, p2] : cell2Sorted._particles) {
        if (std::abs(p1Projection - p2Projection) > _sortingCutoff) {
          break;
        }
        interactParticlesNoN3(*p1, *p2);
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

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCellPairSoAN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3,
                 bidirectional>::processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if constexpr (bidirectional) {
    _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCellSoAN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption::Value dataLayout,
          bool useNewton3, bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, dataLayout, useNewton3, bidirectional>::processCellSoANoN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}
}  // namespace autopas::internal
