/**
 * @file DirectSumContainer.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/inBox.h"

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
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff, allDSApplicableTraversals()),
        _cellBorderFlagManager() {
    this->_cells.resize(2);
  }

  /**
   * Lists all traversal options applicable for the Direct Sum container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOptions> &allDSApplicableTraversals() {
    static const std::vector<TraversalOptions> v{TraversalOptions::directSumTraversal};
    return v;
  }

  ContainerOptions getContainerType() override { return ContainerOptions::directSumContainer; }

  void addParticle(Particle &p) override {
    if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      getCell()->addParticle(p);
    } else {
      utils::ExceptionHandler::exception(
          "DirectSum: trying to add particle that is not in the bounding box.\n" + p.toString());
    }
  }

  void addHaloParticle(Particle &p) override {
    if (utils::notInBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      getHaloCell()->addParticle(p);
    } else {  // particle is not outside of own box
      utils::ExceptionHandler::exception(
          "DirectSum: trying to add particle that is not OUTSIDE of the "
          "bounding box.\n" +
          p.toString());
    }
  }

  void deleteHaloParticles() override { getHaloCell()->clear(); }

  /**
   * @copydoc LinkedCells::iteratePairwiseAoS
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = true) {
    AutoPasLog(debug, "Using traversal {} with AoS", traversal->getTraversalType());
    traversal->traverseCellPairs(this->_cells);
  }

  /**
   * @copydoc LinkedCells::iteratePairwiseSoA
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseSoA(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = true) {
    AutoPasLog(debug, "Using traversal {} with SoA ", traversal->getTraversalType());
    f->SoALoader(*getCell(), (*getCell())._particleSoABuffer);
    f->SoALoader(*getHaloCell(), (*getHaloCell())._particleSoABuffer);

    traversal->traverseCellPairs(this->_cells);

    f->SoAExtractor((*getCell()), (*getCell())._particleSoABuffer);
    f->SoAExtractor((*getHaloCell()), (*getHaloCell())._particleSoABuffer);
  }

  void updateContainer() override {
    if (getHaloCell()->isNotEmpty()) {
      utils::ExceptionHandler::exception(
          "DirectSum: Halo particles still present when updateContainer was called. Found {} particles",
          getHaloCell()->numParticles());
    }
    for (auto iter = getCell()->begin(); iter.isValid(); ++iter) {
      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        addHaloParticle(*iter);
        iter.deleteCurrentParticle();
      }
    }
  }

  bool isContainerUpdateNeeded() override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for ???
#pragma omp parallel shared(outlierFound)  // if (this->_cells.size() / omp_get_max_threads() > ???)
#endif
    for (auto iter = this->begin(); iter.isValid() && (not outlierFound); ++iter) {
      if (utils::notInBox(iter->getR(), this->getBoxMin(), this->getBoxMax())) {
        outlierFound = true;
      }
    }
    return outlierFound;
  }

  TraversalSelector<ParticleCell> generateTraversalSelector(std::vector<TraversalOptions> traversalOptions) override {
    std::vector<TraversalOptions> allowedAndApplicable;

    std::sort(traversalOptions.begin(), traversalOptions.end());
    std::set_intersection(this->_applicableTraversals.begin(), this->_applicableTraversals.end(),
                          traversalOptions.begin(), traversalOptions.end(), std::back_inserter(allowedAndApplicable));
    // direct sum technically consists of two cells (owned + halo)
    return TraversalSelector<ParticleCell>({2, 0, 0}, allowedAndApplicable);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_cells, 0, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(std::array<double, 3> lowerCorner,
                                                      std::array<double, 3> higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    std::vector<size_t> cellsOfInterest;

    switch (behavior) {
      case IteratorBehavior::haloOnly:
        cellsOfInterest.push_back(1);
        break;
      case IteratorBehavior::ownedOnly:
        cellsOfInterest.push_back(0);
        break;
      case IteratorBehavior::haloAndOwned:
        cellsOfInterest.push_back(0);
        cellsOfInterest.push_back(1);
        break;
    }

    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, ParticleCell>(
        &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBorderFlagManager, behavior));
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

  ParticleCell *getCell() { return &(this->_cells.at(0)); };

  ParticleCell *getHaloCell() { return &(this->_cells.at(1)); };
};

}  // namespace autopas