/**
 * @file LinkedCells.h
 *
 * @author tchipevn
 * @date 17.02.2018
 */

#pragma once

#include "CellBlock3D.h"
#include "ParticleContainer.h"
#include "containers/cellPairTraversals/C08Traversal.h"
#include "containers/cellPairTraversals/SlicedTraversal.h"
#include "pairwiseFunctors/CellFunctor.h"
#include "utils/inBox.h"
#include "iterators/ParticleIterator.h"
#include "iterators/RegionParticleIterator.h"

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
 private:
  static const std::vector<TraversalOptions> &LCApplicableTraversals() {
    static const std::vector<TraversalOptions> v{TraversalOptions::c08, TraversalOptions::sliced};
    return v;
  }

 public:
  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param allowedTraversalOptions Traversal options from which the TraversalSelector shall choose.
   * By default all applicable traversals are allowed.
   */
  LinkedCells(const std::array<double, 3> boxMin,
              const std::array<double, 3> boxMax,
              const double cutoff,
              const unsigned int retuneInterval = 100,
              const std::vector<TraversalOptions> &allowedTraversalOptions = LCApplicableTraversals())
      : ParticleContainer<Particle, ParticleCell>(boxMin,
                                                  boxMax,
                                                  cutoff,
                                                  LCApplicableTraversals()),
        _cellBlock(this->_data,
                   boxMin,
                   boxMax,
                   cutoff) {

    assert(this->checkIfTraversalsAreApplicable(allowedTraversalOptions));
    if (not allowedTraversalOptions.empty())
      this->_traversalSelector =
          new TraversalSelector<ParticleCell>(_cellBlock.getCellsPerDimensionWithHalo(),
                                              retuneInterval,
                                              allowedTraversalOptions);

  }

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

  void iteratePairwiseAoS(Functor<Particle, ParticleCell> *f, bool useNewton3 = true) override {
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
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, true> cellFunctor(f);
      this->_traversalSelector->getOptimalTraversal(cellFunctor)->traverseCellPairs(this->_data);
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, false> cellFunctor(f);
      this->_traversalSelector->getOptimalTraversal(cellFunctor)->traverseCellPairs(this->_data);
    }
  }

  void iteratePairwiseSoA(Functor<Particle, ParticleCell> *f, bool useNewton3 = true) override {
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
    loadSoAs(f);

    if (useNewton3) {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, true> cellFunctor(f);
      this->_traversalSelector->getOptimalTraversal(cellFunctor)->traverseCellPairs(this->_data);
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, false> cellFunctor(f);
      this->_traversalSelector->getOptimalTraversal(cellFunctor)->traverseCellPairs(this->_data);
    }

    extractSoAs(f);
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

  bool isContainerUpdateNeeded() override {
    for (size_t cellIndex1d = 0; cellIndex1d < this->_data.size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _cellBlock.getCellBoundingBox(cellIndex1d, boxmin, boxmax);
      for (auto iter = this->_data[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (not inBox(iter->getR(), boxmin, boxmax)) {
          return true;  // we need an update
        }
      }
    }
    return false;
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_data, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, ParticleCell>(
        &this->_data, lowerCorner, higherCorner, &_cellBlock, behavior));
  }

 protected:
  /**
   * object to manage the block of cells.
   */
  CellBlock3D<ParticleCell> _cellBlock;
  // ThreeDimensionalCellHandler

  /**
   * Iterate over all cells and load the data in the SoAs.
   * @tparam ParticleFunctor
   * @param functor
   */
  template <class ParticleFunctor>
  void loadSoAs(ParticleFunctor *functor) {
#ifdef AUTOPAS_OPENMP
    // TODO find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < this->_data.size(); ++i) {
      functor->SoALoader(this->_data[i], this->_data[i]._particleSoABuffer);
    }
  }

  /**
   * Iterate over all cells and fetch the data from the SoAs.
   * @tparam ParticleFunctor
   * @param functor
   */
  template <class ParticleFunctor>
  void extractSoAs(ParticleFunctor *functor) {
#ifdef AUTOPAS_OPENMP
    // TODO find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < this->_data.size(); ++i) {
      functor->SoAExtractor(this->_data[i], this->_data[i]._particleSoABuffer);
    }
  }
};
}  // namespace autopas
