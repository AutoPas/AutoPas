/*
 * Direct.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_DIRECTSUM_H_
#define SRC_CONTAINERS_DIRECTSUM_H_

#include "ParticleContainer.h"
#include "pairwiseFunctors/CellFunctor.h"
#include "pairwiseFunctors/LJFunctor.h"
#include "utils/inBox.h"

namespace autopas {

template <class Particle, class ParticleCell>
class DirectSum : public ParticleContainer<Particle, ParticleCell> {
 public:
  DirectSum(const std::array<double, 3> boxMin,
            const std::array<double, 3> boxMax, double cutoff)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff) {
    this->_data.resize(2);
  }

  void addParticle(Particle &p) override {
    bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      getCell()->addParticle(p);
    } else {
      // todo
    }
  }

  void addHaloParticle(Particle &p) override {
    bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (not inBox) {
      getHaloCell()->addParticle(p);
    } else {  // particle is not outside of own box
      // todo
    }
  }

  void deleteHaloParticles() override {
    getHaloCell()->clear();
  }

  void iteratePairwiseAoS(Functor<Particle, ParticleCell> *f) override {
    //		CellFunctor<Particle, ParticleCell,LJFunctor<Particle>>
    // cellFunctor(f);
    //		cellFunctor.processCellAoSN3(*getCell());
    iteratePairwiseAoS2(f);
#if 0
		for (auto outer = getIt(); outer.isValid(); ++outer) {
			Particle & p1 = *outer;

			int ind = outer.getIndex() + 1;

			for (auto inner = getIt(ind); inner.isValid(); ++inner) {
				Particle & p2 = *inner;

				f->AoSFunctor(p1, p2);
			}
		}
#endif
  }

  template <class ParticleFunctor>
  void iteratePairwiseAoS2(ParticleFunctor *f) {
    CellFunctor<Particle, ParticleCell, ParticleFunctor, false> cellFunctor(f);
    cellFunctor.processCell(*getCell());
  }

  void iteratePairwiseSoA(Functor<Particle, ParticleCell> *f) override {
    iteratePairwiseSoA2(f);
  }

  template <class ParticleFunctor>
  void iteratePairwiseSoA2(ParticleFunctor *f) {
    CellFunctor<Particle, ParticleCell, ParticleFunctor, true> cellFunctor(f);
    cellFunctor.processCell(*getCell());
  }

  void updateContainer() override {
    // TODO: might need to do sth. if particles move outside of the box?
  }

 private:

  // for convenience
  typedef SingleCellIterator<Particle, ParticleCell> SingIterator;

  SingIterator getIt(int index = 0) { return SingIterator(getCell(), index); }

  ParticleCell* getCell() { return &(this->_data.at(0)); };

  ParticleCell* getHaloCell() { return &(this->_data.at(1)); };
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_DIRECTSUM_H_ */
