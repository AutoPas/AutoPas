/**
 * @file DirectSum.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/CudaStreamHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/StringUtils.h"
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
  using ParticleFloatType = typename Particle::ParticleFloatingPointType;

 public:
  /**
   * Constructor of the DirectSum class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   */
  DirectSum(const std::array<ParticleFloatType, 3> boxMin, const std::array<ParticleFloatType, 3> boxMax,
            ParticleFloatType cutoff)
      : ParticleContainer<Particle, ParticleCell>(boxMin, boxMax, cutoff, allDSApplicableTraversals()),
        _cellBorderFlagManager() {
    this->_cells.resize(2);
  }

  /**
   * Lists all traversal options applicable for the Direct Sum container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOption> &allDSApplicableTraversals() {
    static const std::vector<TraversalOption> v{TraversalOption::directSumTraversal};
    return v;
  }

  std::vector<TraversalOption> getAllTraversals() override { return allDSApplicableTraversals(); }

  ContainerOption getContainerType() override { return ContainerOption::directSum; }

  void addParticle(Particle &p) override {
    if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      getCell()->addParticle(p);
    } else {
      utils::ExceptionHandler::exception("DirectSum: trying to add a particle that is not in the bounding box.\n" +
                                         p.toString());
    }
  }

  void addHaloParticle(Particle &p) override {
    if (utils::notInBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
      getHaloCell()->addParticle(p);
    } else {  // particle is not outside of own box
      utils::ExceptionHandler::exception(
          "DirectSum: trying to add a halo particle that is not OUTSIDE of the "
          "bounding box.\n" +
          p.toString());
    }
  }

  void deleteHaloParticles() override { getHaloCell()->clear(); }

  /**
   * Function to iterate over all pairs of particles
   * @tparam ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = false) {
    AutoPasLog(debug, "Using traversal {}.", utils::StringUtils::to_string(traversal->getTraversalType()));

    traversal->initTraversal(this->_cells);
    if (auto *traversalInterface = dynamic_cast<DirectSumTraversalInterface<ParticleCell> *>(traversal)) {
      traversalInterface->traverseCellPairs(this->_cells);

    } else {
      autopas::utils::ExceptionHandler::exception(
          "trying to use a traversal of wrong type in DirectSum::iteratePairwise");
    }
    traversal->endTraversal(this->_cells);
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

  TraversalSelector<ParticleCell> generateTraversalSelector() override {
    // direct sum technically consists of two cells (owned + halo)
    return TraversalSelector<ParticleCell>({2, 0, 0});
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_cells, 0, &_cellBorderFlagManager, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(std::array<ParticleFloatType, 3> lowerCorner,
                                                      std::array<ParticleFloatType, 3> higherCorner,
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
