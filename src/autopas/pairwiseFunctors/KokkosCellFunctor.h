/**
 * @file KokkosCellFunctor.h
 *
 * @date 8 Dec 2021
 * @author lgaertner
 */

#pragma once

#include "autopas/cells/KokkosParticleCell.h"
#include "autopas/cells/SortedCellView.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::internal {
/**
 * A cell functor. This functor is build from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles.
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam DataLayout the DataLayout to be used
 * @tparam useNewton3
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3 = true,
          bool bidirectional = true>
class KokkosCellFunctor {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f The ParticleFunctor which should be used for the interaction.
   * @param sortingCutoff This paramater indicates the maximal distance the sorted particles are to interact. This
   * parameter is only relevant for optimization (sorting). This parameter normally should be the cutoff, for building
   * verlet lists, this should be cutoff+skin.
   */
  explicit KokkosCellFunctor(ParticleFunctor *f) : _functor(f) {}

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
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2,
                            const std::array<double, 3> &sortingDirection);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param sortingDirection Normalized vector connecting centers of cell1 and cell2.
   */
  void processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                              const std::array<double, 3> &sortingDirection);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  /**
   * Min. number of particles to start sorting.
   * @todo Currently, this is disabled because of https://github.com/AutoPas/AutoPas/issues/418
   */
  constexpr static unsigned long _startSorting = std::numeric_limits<unsigned long>::max();
};

template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCell(
    ParticleCell &cell) {
  if ((DataLayout == DataLayoutOption::soa && cell._particleSoABuffer.getNumParticles() == 0) ||
      (DataLayout == DataLayoutOption::aos && cell.numParticles() == 0)) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      processCellAoS<useNewton3>(cell);
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
        processCellSoAN3(cell);
      } else {
        processCellSoANoN3(cell);
      }
      break;
  }
}

template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPair(

    ParticleCell &cell1, ParticleCell &cell2,
    const std::array<double, 3> &sortingDirection) {
  if ((DataLayout == DataLayoutOption::soa &&
       (cell1._particleSoABuffer.getNumParticles() == 0 || cell2._particleSoABuffer.getNumParticles() == 0)) ||
      (DataLayout == DataLayoutOption::aos && (cell1.numParticles() == 0 || cell2.numParticles() == 0))) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if (useNewton3) {
        processCellPairAoSN3(cell1, cell2, sortingDirection);
      } else {
        processCellPairAoSNoN3(cell1, cell2, sortingDirection);
      }
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
      break;
  }
}

template <class ParticleCell, class ParticleFunctor, autopas::DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
template <bool newton3>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellAoS(
    ParticleCell &cell) {
  // NOTE: no sorting within cell
  size_t outer = 0ul;
  for (; outer < cell.getSize(); ++outer) {
    auto inner = outer;
    ++inner;
    for (; inner < cell.getSize(); ++inner) {
      if constexpr (newton3) {
        _functor->AoSFunctor(cell[outer], cell[inner], true);
      } else {
        _functor->AoSFunctor(cell[inner], cell[outer], false);
        _functor->AoSFunctor(cell[outer], cell[inner], false);
      }
    }
  }
}

template <class ParticleCell, class ParticleFunctor, autopas::DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2,
    const std::array<double, 3> &sortingDirection) {
  // NOTE: no sorting within cell
  for (size_t outer = 0; outer < cell1.getSize(); outer++) {
    for (size_t inner = 0; inner < cell2.getSize(); inner++) {
      _functor->AoSFunctor(cell1[outer], cell2[inner], true);
    }
  }
}

template <class ParticleCell, class ParticleFunctor, autopas::DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairAoSNoN3(
    ParticleCell &cell1, ParticleCell &cell2,
    const std::array<double, 3> &sortingDirection) {
  // NOTE: no sorting within cell
  for (size_t outer = 0; outer < cell1.getSize(); outer++) {
    for (size_t inner = 0; inner < cell2.getSize(); inner++) {
      _functor->AoSFunctor(cell1[outer], cell2[inner], false);
      if (bidirectional) _functor->AoSFunctor(cell2[inner], cell1[outer], false);
    }
  }
}

template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairSoAN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairSoANoN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if (bidirectional) _functor->SoAFunctorPair(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
}

template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoAN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, true);
}

template <class ParticleCell, class ParticleFunctor, DataLayoutOption::Value DataLayout, bool useNewton3,
          bool bidirectional>
void KokkosCellFunctor<ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoANoN3(
    ParticleCell &cell) {
  _functor->SoAFunctorSingle(cell._particleSoABuffer, false);  // the functor has to enable this...
}
}  // namespace autopas::internal
