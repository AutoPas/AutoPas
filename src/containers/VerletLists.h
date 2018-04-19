/**
 * @file VerletLists.h
 * @author seckler
 * @date 19.04.18
 */

#pragma once

#include "LinkedCells.h"

namespace autopas {

/**
 * Verlet Lists container.
 * This class builds neighbour lists for the particle interactions.
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class VerletLists : public LinkedCells<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the VerletLists class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   */
  VerletLists(const std::array<double, 3> boxMin,
              const std::array<double, 3> boxMax, double cutoff)
      : LinkedCells<Particle, ParticleCell>(boxMin, boxMax, cutoff) {}

  void iteratePairwiseAoS(Functor<Particle, ParticleCell> *f,
                          bool useNewton3 = true) override {
    iteratePairwiseAoS2(f, useNewton3);
  }

  /**
   * same as iteratePairwiseAoS, but potentially faster (if called with the
   * derived functor), as the class of the functor is known and thus the
   * compiler can do some better optimizations.
   * @tparam ParticleFunctor
   * @param f
   * @param useNewton3 defines whether newton3 should be used
   */
  template <class ParticleFunctor>
  void iteratePairwiseAoS2(ParticleFunctor *f, bool useNewton3 = true) {
    if (useNewton3) {
      this->updateVerletListsN3();
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, true>
          cellFunctor(f);
      //		cellFunctor.processCellAoSN3(this->_data[13]);
      SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell,
                                                ParticleFunctor, false, true>>
          traversal(this->_data,
                    this->_cellBlock.getCellsPerDimensionWithHalo(),
                    &cellFunctor);

      traversal.traverseCellPairs();
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, false>
          cellFunctor(f);
      //		cellFunctor.processCellAoSN3(this->_data[13]);
      SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell,
                                                ParticleFunctor, false, false>>
          traversal(this->_data,
                    this->_cellBlock.getCellsPerDimensionWithHalo(),
                    &cellFunctor);

      traversal.traverseCellPairs();
    }
  }

 private:
  class VerletListGeneratorFunctor
      : public autopas::Functor<Particle, ParticleCell> {};

  void updateVerletListsN3() {
    VerletListGeneratorFunctor f;
    CellFunctor<Particle, ParticleCell, VerletListGeneratorFunctor, false, true>
        cellFunctor(&f);
    //		cellFunctor.processCellAoSN3(this->_data[13]);
    SlicedTraversal<ParticleCell,
                    CellFunctor<Particle, ParticleCell,
                                VerletListGeneratorFunctor, false, true>>
        traversal(this->_data, this->_cellBlock.getCellsPerDimensionWithHalo(),
                  &cellFunctor);

    traversal.traverseCellPairs();
  }
  // ThreeDimensionalCellHandler
};

} /* namespace autopas */
