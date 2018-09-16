/**
 * @file CellFunctor.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * A cell functor. This functor is build from the normal Functor of the template
 * type ParticleFunctor. It is an internal object to handle interactions between
 * two cells of particles.
 * @todo: currently always used newton3!
 * @tparam Particle
 * @tparam ParticleCell
 * @tparam ParticleFunctor the functor which
 * @tparam useSoA
 * @tparam useNewton3
 * @tparam bidirectional if no newton3 is used processCellPair(cell1, cell2) should also handle processCellPair(cell2,
 * cell1)
 */
template <class Particle, class ParticleCell, class ParticleFunctor, bool useSoA, bool useNewton3 = true,
          bool bidirectional = true>
class CellFunctor {
 public:
  /**
   * The constructor of CellFunctor
   * @param f the particlefunctor which should be used for the interaction.
   */
  explicit CellFunctor(ParticleFunctor *f) : _functor(f) {}

  /**
   * process the interactions inside one cell
   * @param cell all pairwise interactions of particles inside this cell are
   * calculated
   */
  void processCell(ParticleCell &cell) {
    if (cell.numParticles() == 0) {
      return;
    }
    if (useSoA) {
      if (useNewton3) {
        processCellSoAN3(cell);
      } else {
        processCellSoANoN3(cell);
      }
    } else {
      if (useNewton3) {
        processCellAoSN3(cell);
      } else {
        processCellAoSNoN3(cell);
      }
    }
  }

  /**
   * process the interactions between the particles of cell1 with particles of
   * cell2.
   * @param cell1
   * @param cell2
   */
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2) {
    if (cell1.numParticles() == 0 || cell2.numParticles() == 0) {
      return;
    }
    if (useSoA) {
      if (useNewton3) {
        processCellPairSoAN3(cell1, cell2);
      } else {
        processCellPairSoANoN3(cell1, cell2);
      }
    } else {
      if (useNewton3) {
        processCellPairAoSN3(cell1, cell2);
      } else {
        processCellPairAoSNoN3(cell1, cell2);
      }
    }
  }

 private:
  /**
   * Applies the functor to all particle pairs exploiting newtons third law of
   * motion
   * @param cell
   */
  void processCellAoSN3(ParticleCell &cell) {
    AUTOPAS_WITH_STATIC_CELL_ITER(outer, cell, {
      for (; outer.isValid(); ++outer) {
        Particle &p1 = *outer;

        auto inner = outer;
        ++inner;
        for (; inner.isValid(); ++inner) {
          Particle &p2 = *inner;

          _functor->AoSFunctor(p1, p2, true);
        }
      }
    })
  }

  /**
   * Applies the functor to all particle pairs without exploiting newtons third
   * law of motion
   * @param cell
   */
  void processCellAoSNoN3(ParticleCell &cell) {
    AUTOPAS_WITH_STATIC_CELL_ITER(outer, cell, {
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
    })
  }

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion
   * @param cell1
   * @param cell2
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2) {
    AUTOPAS_WITH_STATIC_CELL_ITER(outer, cell1, {
      AUTOPAS_WITH_STATIC_CELL_ITER(innerStart, cell2, {
        // body
        for (; outer.isValid(); ++outer) {
          Particle &p1 = *outer;

          for (auto inner = innerStart; inner.isValid(); ++inner) {
            Particle &p2 = *inner;

            _functor->AoSFunctor(p1, p2, true);
          }
        }
      });
    });
  }

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * without exploiting newtons third law of motion
   * @param cell1
   * @param cell2
   */
  void processCellPairAoSNoN3(ParticleCell &cell1, ParticleCell &cell2) {
    AUTOPAS_WITH_STATIC_CELL_ITER(outer, cell1, {
      AUTOPAS_WITH_STATIC_CELL_ITER(innerStart, cell2, {
        // body
        for (auto outer = cell1.begin(); outer.isValid(); ++outer) {
          Particle &p1 = *outer;

          for (auto inner = innerStart; inner.isValid(); ++inner) {
            Particle &p2 = *inner;

            _functor->AoSFunctor(p1, p2, false);
            if (bidirectional) _functor->AoSFunctor(p2, p1, false);
          }
        }
      });
    });
  }

  void processCellPairSoAN3(ParticleCell &cell1, ParticleCell &cell2) {
    _functor->SoAFunctor(cell1._particleSoABuffer, cell2._particleSoABuffer, true);
  }

  void processCellPairSoANoN3(ParticleCell &cell1, ParticleCell &cell2) {
    _functor->SoAFunctor(cell1._particleSoABuffer, cell2._particleSoABuffer, false);
    if (bidirectional) _functor->SoAFunctor(cell2._particleSoABuffer, cell1._particleSoABuffer, false);
  }

  void processCellSoAN3(ParticleCell &cell) { _functor->SoAFunctor(cell._particleSoABuffer, true); }

  void processCellSoANoN3(ParticleCell &cell) {
    _functor->SoAFunctor(cell._particleSoABuffer,
                         false);  // the functor has to enable this...
  }

 private:
  ParticleFunctor *_functor;
};

}  // namespace autopas
