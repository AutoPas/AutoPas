/**
 * @file CellFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/cells/FullSortedParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {
namespace internal {
/**
 * A cell functor. This functor is build from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles.
 * @todo: currently always used newton3!
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam DataLayout the DataLayout to be used
 * @tparam useNewton3
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout,
          bool useNewton3 = true, bool bidirectional = true>
class CellFunctor {
 public:
  /**
   * The constructor of CellFunctor.
   * @param f the particlefunctor which should be used for the interaction.
   * @param cutoff cutoff radius.
   */
  explicit CellFunctor(ParticleFunctor *f, const double cutoff = std::numeric_limits<double>::max())
      : _functor(f), _cutoff(cutoff) {}

  /**
   * process the interactions inside one cell.
   * @param cell all pairwise interactions of particles inside this cell are
   * calculated
   */
  void processCell(ParticleCell &cell);

  /**
   * process the interactions between the particles of cell1 with particles of
   * cell2.
   * @param cell1
   * @param cell2
   * @param r normalized vector connecting centers of cell1 and cell2
   */
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &r = {1., 0., 0.});

 private:
  /**
   * Applies the functor to all particle pairs exploiting newtons third law of
   * motion.
   * @param cell
   */
  void processCellAoSN3(ParticleCell &cell);

  /**
   * Applies the functor to all particle pairs without exploiting newtons third
   * law of motion.
   * @param cell
   */
  void processCellAoSNoN3(ParticleCell &cell);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param r normalized vector connecting centers of cell1 and cell2
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &r);

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion.
   * @param cell1
   * @param cell2
   * @param r normalized vector connecting centers of cell1 and cell2
   */

  void processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &r);

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairCudaNoN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellPairCudaN3(ParticleCell &cell1, ParticleCell &cell2);

  void processCellSoAN3(ParticleCell &cell);

  void processCellSoANoN3(ParticleCell &cell);

  void processCellCudaNoN3(ParticleCell &cell);

  void processCellCudaN3(ParticleCell &cell);

  ParticleFunctor *_functor;

  const double _cutoff;
};

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCell(
    ParticleCell &cell) {
  if ((DataLayout == DataLayoutOption::soa && cell._particleSoABuffer.getNumParticles() == 0) ||
      (DataLayout == DataLayoutOption::aos && cell.numParticles() == 0)) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if (useNewton3) {
        processCellAoSN3(cell);
      } else {
        processCellAoSNoN3(cell);
      }
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
        processCellSoAN3(cell);
      } else {
        processCellSoANoN3(cell);
      }
      break;
    case DataLayoutOption::cuda:
      if (useNewton3) {
        processCellCudaN3(cell);
      } else {
        processCellCudaNoN3(cell);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPair(

    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &r) {
  if ((DataLayout == DataLayoutOption::soa &&
       (cell1._particleSoABuffer.getNumParticles() == 0 || cell2._particleSoABuffer.getNumParticles() == 0)) ||
      (DataLayout == DataLayoutOption::aos && (cell1.numParticles() == 0 || cell2.numParticles() == 0))) {
    return;
  }

  switch (DataLayout) {
    case DataLayoutOption::aos:
      if (useNewton3) {
        processCellPairAoSN3(cell1, cell2, r);
      } else {
        processCellPairAoSNoN3(cell1, cell2, r);
      }
      break;
    case DataLayoutOption::soa:
      if (useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
      break;
    case DataLayoutOption::cuda:
      if (useNewton3) {
        processCellPairCudaN3(cell1, cell2);
      } else {
        processCellPairCudaNoN3(cell1, cell2);
      }
      break;
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellAoSN3(
    ParticleCell &cell) {
  if (cell.numParticles() > 8) {
    FullSortedParticleCell<Particle, ParticleCell> cellSorted(
        cell, ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    auto outer = cellSorted._particles.begin();
    for (; outer != cellSorted._particles.end(); ++outer) {
      Particle &p1 = *outer->second;

      auto inner = outer;
      ++inner;
      for (; inner != cellSorted._particles.end(); ++inner) {
        if (std::abs(outer->first - inner->first) > _cutoff) {
          break;
        }
        Particle &p2 = *inner->second;
        _functor->AoSFunctor(p1, p2, true);
      }
    }
  } else {
    auto outer = getStaticCellIter(cell);
    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      auto inner = outer;
      ++inner;
      for (; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2, true);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellAoSNoN3(
    ParticleCell &cell) {
  if (cell.numParticles() > 8) {
    FullSortedParticleCell<Particle, ParticleCell> cellSorted(
        cell, ArrayMath::normalize(std::array<double, 3>{1.0, 1.0, 1.0}));

    auto outer = cellSorted._particles.begin();
    auto innerStart = outer;

    for (; outer != cellSorted._particles.end(); ++outer) {
      Particle &p1 = *outer->second;

      // loop over everything until outer
      auto inner = innerStart;
      for (; inner != outer; ++inner) {
        Particle &p2 = *inner->second;

        _functor->AoSFunctor(p1, p2, false);
      }
      // skip over the outer one
      ++inner;

      // loop over everything after outer
      for (; inner != cellSorted._particles.end(); ++inner) {
        if (std::abs(outer->first - inner->first) > _cutoff) {
          break;
        }
        Particle &p2 = *inner->second;
        _functor->AoSFunctor(p1, p2, false);
      }
    }
  } else {
    auto outer = getStaticCellIter(cell);
    auto innerStart = outer;

    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      // loop over everything until outer
      auto inner = innerStart;
      for (; inner != outer; ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2, false);
      }
      // skip over the outer one
      ++inner;

      // loop over everything after outer
      for (; inner.isValid(); ++inner) {
        Particle &p2 = *inner;
        _functor->AoSFunctor(p1, p2, false);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairAoSN3(
    ParticleCell &cell1, ParticleCell &cell2, const std::array<double, 3> &r) {
  if (cell1.numParticles() + cell2.numParticles() > 8) {
    FullSortedParticleCell<Particle, ParticleCell> baseSorted(cell1, r);
    FullSortedParticleCell<Particle, ParticleCell> outerSorted(cell2, r);

    for (auto &outer : baseSorted._particles) {
      Particle &p1 = *outer.second;

      for (auto &inner : outerSorted._particles) {
        if (std::abs(outer.first - inner.first) > _cutoff) {
          break;
        }
        Particle &p2 = *inner.second;
        _functor->AoSFunctor(p1, p2, true);
      }
    }
  } else {
    auto outer = getStaticCellIter(cell1);
    auto innerStart = getStaticCellIter(cell2);

    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      for (auto inner = innerStart; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2, true);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                 bidirectional>::processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2,
                                                        const std::array<double, 3> &r) {
  if (cell1.numParticles() + cell2.numParticles() > 8) {
    FullSortedParticleCell<Particle, ParticleCell> baseSorted(cell1, r);
    FullSortedParticleCell<Particle, ParticleCell> outerSorted(cell2, r);

    for (auto &outer : baseSorted._particles) {
      Particle &p1 = *outer.second;

      for (auto &inner : outerSorted._particles) {
        if (std::abs(outer.first - inner.first) > _cutoff) {
          break;
        }
        Particle &p2 = *inner.second;
        _functor->AoSFunctor(p1, p2, false);
        if (bidirectional) _functor->AoSFunctor(p2, p1, false);
      }
    }
  } else {
    auto innerStart = getStaticCellIter(cell2);

    for (auto outer = cell1.begin(); outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      for (auto inner = innerStart; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2, false);
        if (bidirectional) _functor->AoSFunctor(p2, p1, false);
      }
    }
  }
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairSoAN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctor(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                 bidirectional>::processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->SoAFunctor(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
  if (bidirectional) _functor->SoAFunctor(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoAN3(
    ParticleCell &cell) {
  _functor->SoAFunctor(cell._particleSoABuffer, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellSoANoN3(
    ParticleCell &cell) {
  _functor->SoAFunctor(cell._particleSoABuffer, false);  // the functor has to enable this...
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3,
                 bidirectional>::processCellPairCudaNoN3(ParticleCell &cell1, ParticleCell &cell2) {
  _functor->CudaFunctor(cell1._particleSoABufferDevice, cell2._particleSoABufferDevice, false);
  if (bidirectional) _functor->CudaFunctor(cell2._particleSoABufferDevice, cell1._particleSoABufferDevice, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellCudaNoN3(
    ParticleCell &cell) {
  _functor->CudaFunctor(cell._particleSoABufferDevice, false);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellPairCudaN3(
    ParticleCell &cell1, ParticleCell &cell2) {
  _functor->CudaFunctor(cell1._particleSoABufferDevice, cell2._particleSoABufferDevice, true);
}

template <class Particle, class ParticleCell, class ParticleFunctor, DataLayoutOption DataLayout, bool useNewton3,
          bool bidirectional>
void CellFunctor<Particle, ParticleCell, ParticleFunctor, DataLayout, useNewton3, bidirectional>::processCellCudaN3(
    ParticleCell &cell) {
  _functor->CudaFunctor(cell._particleSoABufferDevice, true);
}
}  // namespace internal
}  // namespace autopas
