/*
 * LinkedCells.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef AUTOPAS_SRC_CONTAINERS_LINKEDCELLS_H_
#define AUTOPAS_SRC_CONTAINERS_LINKEDCELLS_H_

#include "CellBlock3D.h"
#include "ParticleContainer.h"
#include "containers/cellPairTraversals/SlicedTraversal.h"
#include "pairwiseFunctors/CellFunctor.h"
#include "utils/inBox.h"

namespace autopas {

/**
 * LinkedCells class.
 * This class uses a list of neighboring cells to store the particles.
 * These cells dimensions at least as large as the given cutoff radius,
 * therefore short-range interactions only need to be calculated between
 * particles in neighboring cells.
 * @tparam Particle type of the particles that need to be stored
 * @tparam ParticleCell type of the ParticleCells that are used to store the
 * particles
 */
template <class Particle, class ParticleCell>
class LinkedCells : public ParticleContainer<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   */
  LinkedCells(const std::array<double, 3> boxMin,
              const std::array<double, 3> boxMax, double cutoff)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff),
        _cellBlock(this->_data, boxMin, boxMax, cutoff) {}

  void addParticle(Particle &p) override {
    bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      ParticleCell &cell = _cellBlock.getContainingCell(p.getR());
      cell.addParticle(p);
    } else {
      utils::ExceptionHandler::exception(
          "LinkedCells: trying to add particle that is not inside the bounding "
          "box");
    }
  }

  void addHaloParticle(Particle &haloParticle) override {
    bool inHalo = _cellBlock.checkInHalo(haloParticle.getR());
    if (inHalo) {
      ParticleCell &cell = _cellBlock.getContainingCell(haloParticle.getR());
      cell.addParticle(haloParticle);
    } else {
      utils::ExceptionHandler::exception(
          "LinkedCells: trying to add halo particle that is not in the halo "
          "box");
    }
  }

  void deleteHaloParticles() override { _cellBlock.clearHaloCells(); }

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
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, true>
          cellFunctor(f);
      //		cellFunctor.processCellAoSN3(this->_data[13]);
      SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell,
                                                ParticleFunctor, false, true>>
          traversal(this->_data, _cellBlock.getCellsPerDimensionWithHalo(),
                    &cellFunctor);

      traversal.traverseCellPairs();
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, false>
          cellFunctor(f);
      //		cellFunctor.processCellAoSN3(this->_data[13]);
      SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell,
                                                ParticleFunctor, false, false>>
          traversal(this->_data, _cellBlock.getCellsPerDimensionWithHalo(),
                    &cellFunctor);

      traversal.traverseCellPairs();
    }
  }

  void iteratePairwiseSoA(Functor<Particle, ParticleCell> *f,
                          bool useNewton3 = true) override {
    /// @todo iteratePairwiseSoA
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
    if (useNewton3) {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, true>
          cellFunctor(f);
      //		cellFunctor.processCellAoSN3(this->_data[13]);
      SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell,
                                                ParticleFunctor, true, true>>
          traversal(this->_data, _cellBlock.getCellsPerDimensionWithHalo(),
                    &cellFunctor);
      traversal.traverseCellPairs();
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, false>
          cellFunctor(f);
      //		cellFunctor.processCellAoSN3(this->_data[13]);
      SlicedTraversal<ParticleCell, CellFunctor<Particle, ParticleCell,
                                                ParticleFunctor, true, false>>
          traversal(this->_data, _cellBlock.getCellsPerDimensionWithHalo(),
                    &cellFunctor);
      traversal.traverseCellPairs();
    }
  }

  void updateContainer() override {
    /// @todo optimize
    std::vector<Particle> invalidParticles;
    for (auto iter = this->begin(); iter.isValid(); ++iter) {
      invalidParticles.push_back(*iter);
    }
    for (auto &cell : this->_data) {
      cell.clear();
    }
    for (auto &particle : invalidParticles) {
      if (inBox(particle.getR(), this->getBoxMin(), this->getBoxMax())) {
        addParticle(particle);
      } else {
        addHaloParticle(particle);
      }
    }
  }

  bool checkUpdateContainerNeeded() override {
    for (int cellIndex1d = 0; cellIndex1d < this->_data.size(); ++cellIndex1d) {
      std::array<double, 3> boxmin;
      std::array<double, 3> boxmax;
      _cellBlock.getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      for (auto iter = this->_data[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (not inBox(iter->getR(),boxmin, boxmax)) {
          return true;  // we need an update
        }
      }
    }
    return false;
  }

 protected:
  /**
   * object to manage the block of cells
   */
  CellBlock3D<ParticleCell> _cellBlock;
  // ThreeDimensionalCellHandler
};

} /* namespace autopas */

#endif /* AUTOPAS_SRC_CONTAINERS_LINKEDCELLS_H_ */
