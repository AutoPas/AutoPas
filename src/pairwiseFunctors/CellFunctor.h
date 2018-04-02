/*
 * CellFunctor.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_

#include "iterators/SingleCellIterator.h"

namespace autopas {

template <class Particle, class ParticleCell, class ParticleFunctor,
          bool useSoA>
class CellFunctor {
 public:
  explicit CellFunctor(ParticleFunctor *f) : _functor(f) {}

  void processCell(ParticleCell &cell) {
    if (useSoA) {
      processCellSoA(cell);
    } else {
      processCellAoSN3(cell);
    }
  }

  void processCellPair(ParticleCell &cell1, ParticleCell &cell2) {
    if (useSoA) {
      processCellPairSoA(cell1, cell2);
    } else {
      processCellPairAoSN3(cell1, cell2);
    }
  }

  /**
   * Applies the functor to all particle pairs exploiting newtons third law of
   * motion
   * @param cell
   */
  void processCellAoSN3(ParticleCell &cell) {
    SingleCellIterator<Particle, ParticleCell> outer(&cell);
    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      int ind = outer.getIndex() + 1;

      SingleCellIterator<Particle, ParticleCell> inner(&cell, ind);
      for (; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2);
      }
    }
  }

  /**
   * Applies the functor to all particle pairs between cell1 and cell2
   * exploiting newtons third law of motion
   * @param cell1
   * @param cell2
   */
  void processCellPairAoSN3(ParticleCell &cell1, ParticleCell &cell2) {
    SingleCellIterator<Particle, ParticleCell> outer(&cell1);
    for (; outer.isValid(); ++outer) {
      Particle &p1 = *outer;

      SingleCellIterator<Particle, ParticleCell> inner(&cell2);
      for (; inner.isValid(); ++inner) {
        Particle &p2 = *inner;

        _functor->AoSFunctor(p1, p2);
      }
    }
  }

  void processCellPairSoA(ParticleCell &cell1, ParticleCell &cell2) {
    _functor->SoALoader(cell1, &cell1._molsSoABuffer);
    _functor->SoALoader(cell2, &cell2._molsSoABuffer);

    _functor->SoAFunctor(cell1._molsSoABuffer, cell2._molsSoABuffer);

    _functor->SoAExtractor(&cell1, &cell1._molsSoABuffer);
    _functor->SoAExtractor(&cell2, &cell2._molsSoABuffer);
  }

  void processCellSoA(ParticleCell &cell) {
    _functor->SoALoader(cell, &cell._molsSoABuffer);

    _functor->SoAFunctor(cell._molsSoABuffer);

    _functor->SoAExtractor(&cell, &cell._molsSoABuffer);
  }

 private:
  ParticleFunctor *_functor;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_ */
