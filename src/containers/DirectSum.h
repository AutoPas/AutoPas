/*
 * Direct.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_DIRECTSUM_H_
#define SRC_CONTAINERS_DIRECTSUM_H_

#include "CellBorderAndFlagManager.h"
#include "ParticleContainer.h"
#include "pairwiseFunctors/CellFunctor.h"
#include "pairwiseFunctors/LJFunctor.h"
#include "utils/ExceptionHandler.h"
#include "utils/inBox.h"

namespace autopas {

/**
 * This class stores particles in a single cell.
 * Interactions are calculated directly, such that each particle interacts with
 * every other particle.
 * Use this class only if you have a very small amount of particles at hand.
 * @tparam Particle type of the particles to be stored
 * @tparam ParticleCell type of the cell that stores the particle
 */
template <class Particle, class ParticleCell>
class DirectSum : public ParticleContainer<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   */
  DirectSum(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, double cutoff)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff) {
    this->_data.resize(2);
  }

  void addParticle(Particle &p) override {
    bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      getCell()->addParticle(p);
    } else {
      utils::ExceptionHandler::exception("DirectSum: trying to add particle that is not in the bounding box");
    }
  }

  void addHaloParticle(Particle &p) override {
    bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (not inBox) {
      getHaloCell()->addParticle(p);
    } else {  // particle is not outside of own box
      utils::ExceptionHandler::exception(
          "DirectSum: trying to add particle that is not OUTSIDE of the "
          "bounding box");
    }
  }

  void deleteHaloParticles() override { getHaloCell()->clear(); }

  void iteratePairwiseAoS(Functor<Particle, ParticleCell> *f, bool useNewton3 = true) override {
    //		CellFunctor<Particle, ParticleCell,LJFunctor<Particle>>
    // cellFunctor(f);
    //		cellFunctor.processCellAoSN3(*getCell());
    iteratePairwiseAoS2(f, useNewton3);
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

  /**
   * same as iteratePairwiseAoS, but faster, as the class of the functor is
   * known and thus the compiler can do some better optimizations.
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3 defines whether newton3 should be used
   */
  template <class ParticleFunctor>
  void iteratePairwiseAoS2(ParticleFunctor *f, bool useNewton3 = true) {
    if (useNewton3) {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, true> cellFunctor(f);
      cellFunctor.processCell(*getCell());
      cellFunctor.processCellPair(*getCell(), *getHaloCell());
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, false> cellFunctor(f);
      cellFunctor.processCell(*getCell());
      cellFunctor.processCellPair(*getCell(), *getHaloCell());
    }
  }

  void iteratePairwiseSoA(Functor<Particle, ParticleCell> *f, bool useNewton3 = true) override {
    iteratePairwiseSoA2(f, useNewton3);
  }

  /**
   * same as iteratePairwiseSoA, but faster, as the class of the functor is
   * known and thus the compiler can do some better optimizations.
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3
   */
  template <class ParticleFunctor>
  void iteratePairwiseSoA2(ParticleFunctor *f, bool useNewton3 = true) {
    f->SoALoader(*getCell(), (*getCell())._particleSoABuffer);
    f->SoALoader(*getHaloCell(), (*getHaloCell())._particleSoABuffer);

    if (useNewton3) {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, true> cellFunctor(f);
      cellFunctor.processCell(*getCell());
      cellFunctor.processCellPair(*getCell(), *getHaloCell());
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, false> cellFunctor(f);
      cellFunctor.processCell(*getCell());
      cellFunctor.processCellPair(*getCell(), *getHaloCell());
    }

    f->SoAExtractor((*getCell()), (*getCell())._particleSoABuffer);
    f->SoAExtractor((*getHaloCell()), (*getHaloCell())._particleSoABuffer);
  }

  void updateContainer() override {
    /// @todo might need to do sth. if particles move outside of the box?
  }

  bool isContainerUpdateNeeded() override {
    for (auto iter = this->begin(); iter.isValid(); ++iter) {
      if (not iter->inBox(this->getBoxMin(), this->getBoxMax())) {
        return true;
      }
    }
    return false;
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_data, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, ParticleCell>(
        &this->_data, lowerCorner, higherCorner, &_cellBorderFlagManager, behavior));
  }

 private:
  class DirectSumCellBorderAndFlagManager : public CellBorderAndFlagManager {
    /**
     * the index type to access the particle cells
     */
    typedef std::size_t index_t;

   public:
    bool isHaloCell(index_t index1d) const override { return index1d == 1; }

    bool isOwningCell(index_t index1d) const override { return not isHaloCell(index1d); }
  } _cellBorderFlagManager;

  ParticleCell *getCell() { return &(this->_data.at(0)); };

  ParticleCell *getHaloCell() { return &(this->_data.at(1)); };
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_DIRECTSUM_H_ */
