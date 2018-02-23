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

template <class Particle, class ParticleCell, class ParticleFunctor>
class CellFunctor {
 public:
  explicit CellFunctor(ParticleFunctor *f) : _functor(f) {}
  void processCell(ParticleCell &cell) { processCellAoSN3(cell); }
  void processCellPair(ParticleCell &cell1, ParticleCell &cell2) {
    processCellPairAoSN3(cell1, cell2);
//    processCellPairSoA(cell1, cell2);
  }

  /**
   * Applies the functor to all particle pairs exploiting newtons third law of motion
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
   * Applies the functor to all particle pairs between cell1 and cell2 exploiting newtons third law of motion
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
    SoA soa1, soa2;

    //TODO: fix vector size (->array)
    std::vector<Particle> particles1, particles2;

    SingleCellIterator<Particle, ParticleCell> iterGet1(&cell1);
    for (; iterGet1.isValid(); ++iterGet1) {
      particles1.push_back(*iterGet1);
    }
    // TODO: Loader sollte aus Zellen laden können
    _functor->SoALoader(particles1, &soa1);

    SingleCellIterator<Particle, ParticleCell> iterGet2(&cell2);
    for (; iterGet2.isValid(); ++iterGet2) {
      particles2.push_back(*iterGet2);
    }
    _functor->SoALoader(particles2, &soa2);

    _functor->SoAFunctor(soa1, soa2);

    particles1.clear();
    // TODO: Extracor sollte in Zellen schreiben können
    _functor->SoAExtractor(&particles1, &soa1);
    particles2.clear();
    _functor->SoAExtractor(&particles2, &soa2);

    SingleCellIterator<Particle, ParticleCell> iterSet1(&cell1);
    for (int i = 0; iterSet1.isValid(); ++iterSet1, ++i) {
      iterSet1->setF(particles1[i].getF());
    }

    SingleCellIterator<Particle, ParticleCell> iterSet2(&cell2);
    for (int i = 0; iterSet2.isValid(); ++iterSet2, ++i) {
      iterSet2->setF(particles2[i].getF());
    }
  }

  void processCellSoA(ParticleCell &cell) {
    SoA soa;

    std::vector<Particle> particles;

    SingleCellIterator<Particle, ParticleCell> iterGet(&cell);
    for (; iterGet.isValid(); ++iterGet) {
      particles.push_back(*iterGet);
    }
    _functor->SoALoader(particles, &soa);

    _functor->SoAFunctor(soa);

    particles.clear();
    _functor->SoAExtractor(&particles, &soa);

    SingleCellIterator<Particle, ParticleCell> iterSet1(&cell);
    for (int i = 0; iterSet1.isValid(); ++iterSet1, ++i) {
      iterSet1->setF(particles[i].getF());
    }
  }

 private:
  ParticleFunctor *_functor;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_ */
