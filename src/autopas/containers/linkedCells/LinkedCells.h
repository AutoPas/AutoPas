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
#include "autopas/particles/FmmParticle.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
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
  using ParticleType = typename ParticleContainer<ParticleCell>::ParticleType;

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

  ContainerOption getContainerType() const override { return ContainerOption::linkedCells; }

  /**
   * @copydoc ParticleContainerInterface::addParticle()
   */
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

  /**
   * @copydoc ParticleContainerInterface::addHaloParticle()
   */
  void addHaloParticle(ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwned(false);
    ParticleCell &cell = _cellBlock.getContainingCell(pCopy.getR());
    cell.addParticle(pCopy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
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
    AutoPasLog(debug, "Using traversal {}.", traversal->getTraversalType().to_string());

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
            internal::deleteParticle(pIter);
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

  bool isContainerUpdateNeeded() const override {
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

  TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(this->getCellBlock().getCellsPerDimensionWithHalo(), this->getInteractionLength(),
                                 this->getCellBlock().getCellLength());
  }

  ParticleIteratorWrapper<ParticleType, true> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::ParticleIterator<ParticleType, ParticleCell, true>(&this->_cells, 0, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<ParticleType, false> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::ParticleIterator<ParticleType, ParticleCell, false>(&this->_cells, 0, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<ParticleType, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::subScalar(lowerCorner, this->getSkin()));
    auto stopIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::addScalar(higherCorner, this->getSkin()));

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

    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, true>(&this->_cells, lowerCorner, higherCorner,
                                                                               cellsOfInterest, &_cellBlock, behavior));
  }

  ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const override {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::subScalar(lowerCorner, this->getSkin()));
    auto stopIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::addScalar(higherCorner, this->getSkin()));

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

    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::RegionParticleIterator<ParticleType, ParticleCell, false>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior));
  }

  /**
   * Get the cell block, not supposed to be used except by verlet lists
   * @return the cell block
   */
  internal::CellBlock3D<ParticleCell> &getCellBlock() { return _cellBlock; }

  /**
   * @copydoc getCellBlock()
   * @note const version
   */
  const internal::CellBlock3D<ParticleCell> &getCellBlock() const { return _cellBlock; }

  /**
   * Returns reference to the data of LinkedCells
   * @return the data
   */
  std::vector<ParticleCell> &getCells() { return this->_cells; }

  /**
   * @copydoc getCells()
   * @note const version
   */
  const std::vector<ParticleCell> &getCells() const { return this->_cells; }

  void createFmmNode(fmm::FmmTreeNode &node, int depth) const {
    // std::cout << "getCellLength " << autopas::utils::ArrayUtils::to_string(_cellBlock.getCellLength()) << std::endl;

    auto halfCell = autopas::utils::ArrayMath::mulScalar(_cellBlock.getCellLength(), 0.5);

    // BoxMin and BoxMax are always between two cells. Adding half a cell moves BoxMin and BoxMax to the center of a
    // cell. This ensures, that the correct 3DIndex is returned.
    auto nodeMin3DIndex = _cellBlock.get3DIndexOfPosition(autopas::utils::ArrayMath::add(node.getBoxMin(), halfCell));
    auto nodeMax3DIndex = _cellBlock.get3DIndexOfPosition(autopas::utils::ArrayMath::add(node.getBoxMax(), halfCell));
    auto delta = utils::ArrayMath::sub(nodeMax3DIndex, nodeMin3DIndex);

    auto tmp = autopas::utils::ArrayMath::sub(nodeMax3DIndex, std::array<unsigned long, 3>({1, 1, 1}));

    node.name = "<<[" + autopas::utils::ArrayUtils::to_string(nodeMin3DIndex) + "] - [" +
                autopas::utils::ArrayUtils::to_string(tmp) + "]>>";

    /*
      std::cout << "-" << std::endl;
      std::cout << "boxMin = " << autopas::utils::ArrayUtils::to_string(node.getBoxMin()) << std::endl;
      std::cout << "boxMax = " << autopas::utils::ArrayUtils::to_string(node.getBoxMax()) << std::endl;
      std::cout << "sphereRadius = " << node.getSphereRadius() << std::endl;

      std::cout << autopas::utils::ArrayUtils::to_string(nodeMin3DIndex) << std::endl;
      std::cout << autopas::utils::ArrayUtils::to_string(nodeMax3DIndex) << std::endl;*/

    // Find the axis in which the box is the largest (based on 3D index).
    int largestIndex = 0;
    if (delta[1] > delta[largestIndex]) {
      largestIndex = 1;
    }
    if (delta[2] > delta[largestIndex]) {
      largestIndex = 2;
    }
    auto largestSize = delta[largestIndex];

    /*std::cout << "delta = " << autopas::utils::ArrayUtils::to_string(delta) << std::endl;
    std::cout << "largestIndex = " << largestIndex << std::endl;
    std::cout << "largestSize = " << largestSize << std::endl;*/

    // Must be at least 2 cells to do a split.
    if (largestSize >= 2) {
      // The 3D index where the node will be split. Only the largestIndex axis is interesting.
      // The other indices are set to the min corner of the box.

      // Split largest size at the center.
      auto largestSizeMid = nodeMin3DIndex[largestIndex] + (largestSize / 2);
      // std::cout << "largestSizeMid = " << largestSizeMid << std::endl;
      auto split3DIndex = nodeMin3DIndex;
      split3DIndex[largestIndex] = largestSizeMid;

      // std::cout << "split3DIndex = " << autopas::utils::ArrayUtils::to_string(split3DIndex) << std::endl;

      // Find the double coordinates of the split position.
      auto midBoxMin = std::array<double, 3>({0, 0, 0});
      auto midBoxMax = std::array<double, 3>({0, 0, 0});
      _cellBlock.getCellBoundingBox(split3DIndex, midBoxMin, midBoxMax);

      // std::cout << "midBoxMin = " << autopas::utils::ArrayUtils::to_string(midBoxMin) << std::endl;
      // std::cout << "midBoxMax = " << autopas::utils::ArrayUtils::to_string(midBoxMax) << std::endl;

      // std::cout << "midBox" << std::endl;

      // std::cout << autopas::ArrayUtils::to_string(midBoxMin) << std::endl;
      // std::cout << autopas::ArrayUtils::to_string(midBoxMax) << std::endl;

      auto splitMin = node.getBoxMin();
      auto splitMax = node.getBoxMax();

      // The lower corner of the center cell is the split position.
      splitMin[largestIndex] = midBoxMin[largestIndex];
      splitMax[largestIndex] = midBoxMin[largestIndex];

      // std::cout << "splitMin = " << autopas::utils::ArrayUtils::to_string(splitMin) << std::endl;
      // std::cout << "splitMax = " << autopas::utils::ArrayUtils::to_string(splitMax) << std::endl;

      // std::cout << autopas::ArrayUtils::to_string(splitMin) << std::endl;
      // std::cout << autopas::ArrayUtils::to_string(splitMax) << std::endl;

      node.split(splitMax, splitMin);
      createFmmNode(node.getChild(0), depth + 1);
      createFmmNode(node.getChild(1), depth + 1);

    } else {
      node.makeLeaf();
    }
  }

  [[nodiscard]] std::unique_ptr<fmm::FmmTree> getFastMultipoleMethodTree() override {
    auto tree = std::make_unique<fmm::FmmTree>();

    createFmmNode(tree->setRoot(this->getBoxMin(), this->getBoxMax()), 0);
    // tree->getRoot().init(tree->getRoot());
    tree->getRoot().initLists(tree->getRoot());
    tree->getRoot().printLists();
    return tree;
  }

 protected:
  /**
   * object to manage the block of cells.
   */
  internal::CellBlock3D<ParticleCell> _cellBlock;
  // ThreeDimensionalCellHandler
};

}  // namespace autopas
