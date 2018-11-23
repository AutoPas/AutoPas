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

#include <bitset>

namespace autopas {

/**
 * LinkedCells class.
 * This class uses a list of neighboring cells to store the particles.
 * These cells dimensions are at least as large as the given cutoff radius,
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
    bool inBox = autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
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
    AutoPasLog(debug, "Using traversal {} with AoS", traversal->getTraversalType());
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
    AutoPasLog(debug, "Using traversal {} with SoA ", traversal->getTraversalType());
    loadSoAs(f);

    traversal->traverseCellPairs(this->_cells);

    extractSoAs(f);
  }

  void updateContainer() override {
    auto haloIter = this->begin(IteratorBehavior::haloOnly);
    if (haloIter.isValid()) {
      utils::ExceptionHandler::exception(
          "Linked Cells: Halo particles still present when updateContainer was called. First particle found:\n" +
          haloIter->toString());
    }

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif  // AUTOPAS_OPENMP
    {
      std::vector<Particle> myInvalidParticles;
#ifdef AUTOPAS_OPENMP
#pragma omp for
#endif  // AUTOPAS_OPENMP
      for (size_t cellId = 0; cellId < this->getCells().size(); ++cellId) {
        // if empty
        if (not this->getCells()[cellId].isNotEmpty()) continue;

        std::array<double, 3> cellLowerCorner, cellUpperCorner;
        this->getCellBlock().getCellBoundingBox(cellId, cellLowerCorner, cellUpperCorner);

        for (auto &&pIter = this->getCells()[cellId].begin(); pIter.isValid(); ++pIter) {
          // if not in cell
          if (utils::notInBox(pIter->getR(), cellLowerCorner, cellUpperCorner)) {
            myInvalidParticles.push_back(*pIter);
            pIter.deleteCurrentParticle();
          }
        }
      }
      // implicit barrier here
      // the barrier is needed because iterators are not threadsafe w.r.t. addParticle()

      // this loop is executed for every thread
      for (auto &&p : myInvalidParticles)
        // if not in halo
        if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax()))
          addParticle(p);
        else
          addHaloParticle(p);
    }
  }

  bool isContainerUpdateNeeded() override {
    std::atomic<bool> outlierFound(false);
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for magic number
    // numThreads should be at least 1 and maximal max_threads
    int numThreads = std::max(1, std::min(omp_get_max_threads(), (int)(this->_cells.size() / 500)));
    AutoPasLog(trace, "Using {} threads", numThreads);
#pragma omp parallel for shared(outlierFound) num_threads(numThreads)
#endif
    for (size_t cellIndex1d = 0; cellIndex1d < this->_cells.size(); ++cellIndex1d) {
      std::array<double, 3> boxmin{0., 0., 0.};
      std::array<double, 3> boxmax{0., 0., 0.};
      _cellBlock.getCellBoundingBox(cellIndex1d, boxmin, boxmax);

      for (auto iter = this->_cells[cellIndex1d].begin(); iter.isValid(); ++iter) {
        if (not utils::inBox(iter->getR(), boxmin, boxmax)) {
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
    std::vector<TraversalOptions> allowedAndApplicable;

    std::sort(traversalOptions.begin(), traversalOptions.end());
    std::set_intersection(this->_applicableTraversals.begin(), this->_applicableTraversals.end(), traversalOptions.begin(),
                          traversalOptions.end(), std::back_inserter(allowedAndApplicable));
    return TraversalSelector<ParticleCell>(this->getCellBlock().getCellsPerDimensionWithHalo(), allowedAndApplicable);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_cells, 0, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(std::array<double, 3> lowerCorner,
                                                      std::array<double, 3> higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    size_t startIndex;
    // this is needed when used through verlet lists since particles can move over cell borders.
    // only lower corner needed since we increase the upper corner anyways.
    if (incSearchRegion) {
      startIndex = this->_cellBlock.get1DIndexOfPosition({
          lowerCorner[0] - 1,
          lowerCorner[1] - 1,
          lowerCorner[2] - 1,
      });
    } else {
      startIndex = this->_cellBlock.get1DIndexOfPosition(lowerCorner);
    }
    auto stopIndex = this->_cellBlock.get1DIndexOfPosition(higherCorner);

    auto startIndex3D =
        utils::ThreeDimensionalMapping::oneToThreeD(startIndex, this->_cellBlock.getCellsPerDimensionWithHalo());
    auto stopIndex3D =
        utils::ThreeDimensionalMapping::oneToThreeD(stopIndex, this->_cellBlock.getCellsPerDimensionWithHalo());

    size_t numCellsOfInterest = (stopIndex3D[0] - startIndex3D[0] + 1) * (stopIndex3D[1] - startIndex3D[1] + 1) *
                                (stopIndex3D[2] - startIndex3D[2] + 1);
    std::vector<size_t> cellsOfInterest(numCellsOfInterest);

    int i = 0;
    for (size_t z = startIndex3D[2]; z <= stopIndex3D[2]; ++z) {
      for (size_t y = startIndex3D[1]; y <= stopIndex3D[1]; ++y) {
        for (size_t x = startIndex3D[0]; x <= stopIndex3D[0]; ++x) {
          cellsOfInterest[i++] =
              utils::ThreeDimensionalMapping::threeToOneD({x, y, z}, this->_cellBlock.getCellsPerDimensionWithHalo());
        }
      }
    }

    return ParticleIteratorWrapper<Particle>(new internal::RegionParticleIterator<Particle, ParticleCell>(
        &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior));
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
