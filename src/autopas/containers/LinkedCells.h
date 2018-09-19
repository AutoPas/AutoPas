/**
 * @file LinkedCells.h
 *
 * @author tchipevn
 * @date 17.02.2018
 */

#pragma once

#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/cellPairTraversals/C08Traversal.h"
#include "autopas/containers/cellPairTraversals/SlicedTraversal.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * LinkedCells class.
 * This class uses a list of neighboring cells to store the particles.
 * These cells dimensions at least as large as the given cutoff radius,
 * therefore short-range interactions only need to be calculated between
 * particles in neighboring cells.
 * @tparam Particle type of the particles that need to be stored
 * @tparam ParticleCell type of the ParticleCells that are used to store the particles
 * @tparam SoAArraysType type of the SoA, needed for verlet lists
 */
template <class Particle, class ParticleCell, class SoAArraysType = typename Particle::SoAArraysType>
class LinkedCells : public ParticleContainer<Particle, ParticleCell, SoAArraysType> {
 public:
  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * By default all applicable traversals are allowed.
   */
  LinkedCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff)
      : ParticleContainer<Particle, ParticleCell, SoAArraysType>(boxMin, boxMax, cutoff, allLCApplicableTraversals()),
        _cellBlock(this->_cells, boxMin, boxMax, cutoff) {}

  /**
   * Lists all traversal options applicable for the Linked Cells container.
   * @return Vector of all applicable traversal options.
   */
  static const std::vector<TraversalOptions> &allLCApplicableTraversals() {
    static const std::vector<TraversalOptions> v{TraversalOptions::c08, TraversalOptions::sliced};
    return v;
  }

  ContainerOptions getContainerType() override { return ContainerOptions::linkedCells; }

  void addParticle(Particle &p) override {
    bool inBox = autopas::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      ParticleCell &cell = _cellBlock.getContainingCell(p.getR());
      cell.addParticle(p);
    } else {
      utils::ExceptionHandler::exception("LinkedCells: trying to add particle that is not inside the bounding box.\n" +
                                         p.toString());
    }
  }

  void addHaloParticle(Particle &haloParticle) override {
    bool inHalo = _cellBlock.checkInHalo(haloParticle.getR());
    if (inHalo) {
      ParticleCell &cell = _cellBlock.getContainingCell(haloParticle.getR());
      cell.addParticle(haloParticle);
    } else {
      utils::ExceptionHandler::exception("LinkedCells: trying to add halo particle that is not in the halo box.\n" +
                                         haloParticle.toString());
    }
  }

  void deleteHaloParticles() override { _cellBlock.clearHaloCells(); }

  /**
   * Function to iterate over all pairs of particles in an array of structures setting. This function only handles
   * short-range interactions.
   * @tparam the type of ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseAoS(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = true) {
    AutoPasLogger->debug("LinkedCells: using traversal {} with AoS", traversal->getTraversalType());
    traversal->traverseCellPairs(this->_cells);
  }

  /**
   * Function to iterate over all pairs of particles in an structure of arrays setting. This function only handles
   * short-range interactions. It is often better vectorizable than iteratePairwiseAoS.
   * @tparam ParticleFunctor
   * @tparam Traversal
   * @param f functor that describes the pair-potential
   * @param traversal the traversal that will be used
   * @param useNewton3 whether newton 3 optimization should be used
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwiseSoA(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = true) {
    AutoPasLogger->debug("LinkedCells: using traversal {} with SoA ", traversal->getTraversalType());
    loadSoAs(f);

    traversal->traverseCellPairs(this->_cells);

    extractSoAs(f);
  }

  void updateContainer() override {
    // @todo optimize
    std::vector<Particle> invalidParticles;
// custom reduction with templates not supported
//#pragma omp parallel reduction(vecMerge: invalidParticles)
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif  // AUTOPAS_OPENMP
    {
      std::vector<Particle> myInvalidParticles;
      for (auto iter = this->begin(); iter.isValid(); ++iter) {
        myInvalidParticles.push_back(*iter);
      }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif  // AUTOPAS_OPENMP
      invalidParticles.insert(invalidParticles.end(), myInvalidParticles.begin(), myInvalidParticles.end());
    }
    for (auto &cell : this->_cells) {
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
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for magic number
    // numThreads should be at least 1 and maximal max_threads
    int numThreads = std::max(1, std::min(omp_get_max_threads(), (int)(this->_cells.size() / 500)));
    AutoPasLogger->trace("LinkedCells::isContainerUpdateNeeded() : using {} threads", numThreads);
#pragma omp parallel for shared(outlierFound) num_threads(numThreads)
#endif
    for (size_t cellIndex1d = 0; cellIndex1d < this->_cells.size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _cellBlock.getCellBoundingBox(cellIndex1d, boxmin, boxmax);

      for (auto iter = this->_cells[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (not inBox(iter->getR(), boxmin, boxmax)) {
          outlierFound = true;  // we need an update
          break;
        }
      }
      // abort loop (for all threads) by moving loop index to end
      if (outlierFound) cellIndex1d = this->_cells.size();
    }

    return outlierFound;
  }

  TraversalSelector<ParticleCell> generateTraversalSelector(std::vector<TraversalOptions> traversalOptions) override {
    return TraversalSelector<ParticleCell>(this->getCellBlock().getCellsPerDimensionWithHalo(), traversalOptions);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_cells, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, ParticleCell>(
        &this->_cells, lowerCorner, higherCorner, &_cellBlock, behavior));
  }

  /**
   * Get the cell block, not supposed to be used except by verlet lists
   * @return the cell block
   */
  CellBlock3D<ParticleCell> &getCellBlock() { return _cellBlock; }

  /**
   * returns reference to the data of LinkedCells
   * @return the data
   */
  std::vector<ParticleCell> &getCells() { return this->_cells; }

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
    // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      functor->SoALoader(this->_cells[i], this->_cells[i]._particleSoABuffer);
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
    // @todo find a condition on when to use omp or when it is just overhead
#pragma omp parallel for
#endif
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      functor->SoAExtractor(this->_cells[i], this->_cells[i]._particleSoABuffer);
    }
  }
};

}  // namespace autopas
