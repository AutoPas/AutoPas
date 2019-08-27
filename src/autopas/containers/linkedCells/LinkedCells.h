/**
 * @file LinkedCells.h
 *
 * @author tchipevn
 * @date 17.02.2018
 */

#pragma once

#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/linkedCells/traversals/LinkedCellTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ParticleCellHelpers.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * LinkedCells class.
 * This class uses a list of neighboring cells to store the particles.
 * These cells dimensions are at least as large as the given cutoff radius,
 * therefore short-range interactions only need to be calculated between
 * particles in neighboring cells.
 * @tparam ParticleCell type of the ParticleCells that are used to store the particles
 * @tparam SoAArraysType type of the SoA, needed for verlet lists
 */
template <class ParticleCell, class SoAArraysType = typename ParticleCell::ParticleType::SoAArraysType>
class LinkedCells : public ParticleContainer<ParticleCell, SoAArraysType> {
 public:
  /**
   *  Type of the Particle.
   */
  typedef typename ParticleCell::ParticleType ParticleType;

  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param cellSizeFactor cell size factor relative to cutoff
   * By default all applicable traversals are allowed.
   */
  LinkedCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
              const double skin, const double cellSizeFactor = 1.0)
      : ParticleContainer<ParticleCell, SoAArraysType>(boxMin, boxMax, cutoff, skin),
        _cellBlock(this->_cells, boxMin, boxMax, cutoff + skin, cellSizeFactor) {}

  ContainerOption getContainerType() override { return ContainerOption::linkedCells; }

  void addParticle(ParticleType &p) override {
    bool inBox = autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      ParticleCell &cell = _cellBlock.getContainingCell(p.getR());
      cell.addParticle(p);
    } else {
      utils::ExceptionHandler::exception(
          "LinkedCells: Trying to add a particle that is not inside the bounding box.\n" + p.toString());
    }
  }

  void addHaloParticle(ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwned(false);
    ParticleCell &cell = _cellBlock.getContainingCell(pCopy.getR());
    cell.addParticle(pCopy);
  }

  bool updateHaloParticle(ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwned(false);
    auto cells = _cellBlock.getNearbyHaloCells(pCopy.getR(), this->getSkin());
    for (auto cellptr : cells) {
      bool updated = internal::checkParticleInCellAndUpdateByID(*cellptr, pCopy);
      if (updated) {
        return true;
      }
    }
    AutoPasLog(trace,
               "UpdateHaloParticle was not able to update particle at "
               "[{}, {}, {}]",
               pCopy.getR()[0], pCopy.getR()[1], pCopy.getR()[2]);
    return false;
  }

  void deleteHaloParticles() override { _cellBlock.clearHaloCells(); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", utils::StringUtils::to_string(traversal->getTraversalType()));

    // Check if traversal is allowed for this container and give it the data it needs.
    auto *traversalInterface = dynamic_cast<LinkedCellTraversalInterface<ParticleCell> *>(traversal);
    auto *cellPairTraversal = dynamic_cast<CellPairTraversal<ParticleCell> *>(traversal);
    if (traversalInterface && cellPairTraversal) {
      cellPairTraversal->setCellsToTraverse(this->_cells);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in LinkedCells::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<ParticleType> updateContainer() override {
    this->deleteHaloParticles();
    std::vector<ParticleType> invalidParticles;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif  // AUTOPAS_OPENMP
    {
      // private for each thread!
      std::vector<ParticleType> myInvalidParticles, myInvalidNotOwnedParticles;
#ifdef AUTOPAS_OPENMP
#pragma omp for
#endif  // AUTOPAS_OPENMP
      for (size_t cellId = 0; cellId < this->getCells().size(); ++cellId) {
        // if empty
        if (not this->getCells()[cellId].isNotEmpty()) continue;

        std::array<double, 3> cellLowerCorner = {}, cellUpperCorner = {};
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

      // this loop is executed for every thread and thus parallel. Don't use #pragma omp for here!
      for (auto &&p : myInvalidParticles) {
        // if not in halo
        if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
          addParticle(p);
        } else {
          myInvalidNotOwnedParticles.push_back(p);
        }
      }
#ifdef AUTOPAS_OPENMP
#pragma omp critical
#endif
      {
        // merge private vectors to global one.
        invalidParticles.insert(invalidParticles.end(), myInvalidNotOwnedParticles.begin(),
                                myInvalidNotOwnedParticles.end());
      }
    }
    return invalidParticles;
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

  TraversalSelectorInfo getTraversalSelectorInfo() override {
    return TraversalSelectorInfo(this->getCellBlock().getCellsPerDimensionWithHalo(), this->getInteractionLength(),
                                 this->getCellBlock().getCellLength());
  }

  ParticleIteratorWrapper<ParticleType> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<ParticleType>(
        new internal::ParticleIterator<ParticleType, ParticleCell>(&this->_cells, 0, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<ParticleType> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D = this->_cellBlock.get3DIndexOfPosition(ArrayMath::subScalar(lowerCorner, this->getSkin()));
    auto stopIndex3D = this->_cellBlock.get3DIndexOfPosition(ArrayMath::addScalar(higherCorner, this->getSkin()));

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

    return ParticleIteratorWrapper<ParticleType>(new internal::RegionParticleIterator<ParticleType, ParticleCell>(
        &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior));
  }

  /**
   * Get the cell block, not supposed to be used except by verlet lists
   * @return the cell block
   */
  internal::CellBlock3D<ParticleCell> &getCellBlock() { return _cellBlock; }

  /**
   * returns reference to the data of LinkedCells
   * @return the data
   */
  std::vector<ParticleCell> &getCells() { return this->_cells; }

 protected:
  /**
   * object to manage the block of cells.
   */
  internal::CellBlock3D<ParticleCell> _cellBlock;
  // ThreeDimensionalCellHandler
};

}  // namespace autopas
