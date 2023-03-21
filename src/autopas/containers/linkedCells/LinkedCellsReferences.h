/**
 * @file LinkedCells.h
 * @date 05.04.2020
 * @author lunaticcoding
 */

#pragma once

#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/cellPairTraversals/BalancedTraversal.h"
#include "autopas/containers/linkedCells/ParticleVector.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/particles/OwnershipState.h"
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
template <class Particle>
class LinkedCellsReferences : public CellBasedParticleContainer<ReferenceParticleCell<Particle>> {
 public:
  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle;
  /**
   *  Type of the ParticleCell.
   */
  using ReferenceCell = ReferenceParticleCell<Particle>;
  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skinPerTimestep
   * @param rebuildFrequency
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   * By default all applicable traversals are allowed.
   */
  LinkedCellsReferences(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                        const double skinPerTimestep, const unsigned int rebuildFrequency,
                        const double cellSizeFactor = 1.0,
                        LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell)
      : CellBasedParticleContainer<ReferenceCell>(boxMin, boxMax, cutoff, skinPerTimestep * rebuildFrequency),
        _cellBlock(this->_cells, boxMin, boxMax, cutoff + skinPerTimestep * rebuildFrequency, cellSizeFactor),
        _loadEstimator(loadEstimator) {std::cout << "actually in linked cell references"<< std::endl;}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::linkedCellsReferences; }

  /**
   * @copydoc ParticleContainerInterface::getParticleCellTypeEnum()
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() override { return CellType::ReferenceParticleCell; }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const ParticleType &p) override {
    ParticleType pCopy = p;
    addParticleLock.lock();
    _particleList.push_back(pCopy);
    updateDirtyParticleReferences();
    addParticleLock.unlock();
  }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);
    addParticleLock.lock();
    _particleList.push_back(pCopy);
    updateDirtyParticleReferences();
    addParticleLock.unlock();
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);
    auto cells = _cellBlock.getNearbyHaloCells(pCopy.getR(), this->getVerletSkin());
    for (auto cellptr : cells) {
      bool updated = internal::checkParticleInCellAndUpdateByID(*cellptr, pCopy);
      if (updated) {
        return true;
      }
    }
    AutoPasLog(trace, "UpdateHaloParticle was not able to update particle: {}", pCopy.toString());
    return false;
  }

  /**
   * @copydoc ParticleContainerInterface::deleteHaloParticles()
   */
  void deleteHaloParticles() override {
    _particleList.clearHaloParticles();
    _cellBlock.clearHaloCells();
  }

  /**
   * Generates the load estimation function depending on _loadEstimator.
   * @return load estimator function object.
   */
  BalancedTraversal::EstimatorFunction getLoadEstimatorFunction() {
    switch (this->_loadEstimator) {
      case LoadEstimatorOption::squaredParticlesPerCell: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          return loadEstimators::squaredParticlesPerCell(this->_cells, cellsPerDimension, lowerCorner, upperCorner);
        };
      }
      case LoadEstimatorOption::none:
        [[fallthrough]];
      default: {
        return
            [&](const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
                const std::array<unsigned long, 3> &upperCorner) { return 1; };
      }
    }
  }

  /**
   * @copydoc ParticleContainerInterface::rebuildNeighborLists()
   */
  void rebuildNeighborLists(TraversalInterface *traversal) override { updateDirtyParticleReferences(); }

  /**
   * Updates all the References in the cells that are out of date.
   * @note Removes all dummy particles.
   */
  void updateDirtyParticleReferences() {
    if (_particleList.isDirty()) {
      if (_particleList.needsRebuild()) {
        for (auto &cell : this->_cells) {
          cell.clear();
        }
      }

      for (auto it = _particleList.beginDirty(); it < _particleList.endDirty(); it++) {
        if (it->isDummy()) {
          continue;
        }
        ReferenceCell &cell = _cellBlock.getContainingCell(it->getR());
        auto address = &(*it);
        cell.addParticleReference(address);
      }
      _particleList.markAsClean();
    }
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    auto *traversalInterface = dynamic_cast<LCTraversalInterface<ReferenceCell> *>(traversal);
    auto *cellPairTraversal = dynamic_cast<CellPairTraversal<ReferenceCell> *>(traversal);
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }
    if (traversalInterface && cellPairTraversal) {
      cellPairTraversal->setCellsToTraverse(this->_cells);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in LinkedCellsReferences::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
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
        // Delete dummy particles of each cell.
        this->getCells()[cellId].deleteDummyParticles();

        // if empty
        if (this->getCells()[cellId].isEmpty()) continue;

        auto [cellLowerCorner, cellUpperCorner] = this->getCellBlock().getCellBoundingBox(cellId);

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
          this->template addParticle<false>(p);
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
    _particleList.deleteDummyParticles();
    updateDirtyParticleReferences();
    return invalidParticles;
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    return TraversalSelectorInfo(this->getCellBlock().getCellsPerDimensionWithHalo(), this->getInteractionLength(),
                                 this->getCellBlock().getCellLength(), 0);
  }

  ParticleIteratorWrapper<ParticleType, true> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) override {
    return ParticleIteratorWrapper<ParticleType, true>(
        new internal::ParticleIterator<ParticleType, ReferenceCell, true>(&this->_cells, 0, &_cellBlock, behavior,
                                                                          nullptr));
  }

  ParticleIteratorWrapper<ParticleType, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const override {
    return ParticleIteratorWrapper<ParticleType, false>(
        new internal::ParticleIterator<ParticleType, ReferenceCell, false>(&this->_cells, 0, &_cellBlock, behavior,
                                                                           nullptr));
  }

  /**
   * @copydoc LinkedCells::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHaloOrDummy) {
    if (behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      // iterate over all particles, so execute directly on particle vector
      _particleList.forEach(forEachLambda);
    } else {
      for (size_t index = 0; index < getCells().size(); index++) {
        if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
          getCells()[index].forEach(forEachLambda, behavior);
        }
      }
    }
  }

  /**
   * @copydoc LinkedCells::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHaloOrDummy) {
    if (behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      // iterate over all particles, so execute directly on particle vector
      _particleList.reduce(reduceLambda, result);
    } else {
      for (size_t index = 0; index < getCells().size(); index++) {
        if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
          getCells()[index].reduce(reduceLambda, result, behavior);
        }
      }
    }
  }

  ParticleIteratorWrapper<ParticleType, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                const std::array<double, 3> &higherCorner,
                                                                IteratorBehavior behavior) override {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::subScalar(lowerCorner, this->getVerletSkin()));
    auto stopIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::addScalar(higherCorner, this->getVerletSkin()));

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
        new internal::RegionParticleIterator<ParticleType, ReferenceCell, true>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior, nullptr));
  }

  ParticleIteratorWrapper<ParticleType, false> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                 const std::array<double, 3> &higherCorner,
                                                                 IteratorBehavior behavior) const override {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::subScalar(lowerCorner, this->getVerletSkin()));
    auto stopIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::addScalar(higherCorner, this->getVerletSkin()));

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
        new internal::RegionParticleIterator<ParticleType, ReferenceCell, false>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior, nullptr));
  }

  /**
   * @copydoc LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::subScalar(lowerCorner, this->getVerletSkin()));
    auto stopIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::addScalar(higherCorner, this->getVerletSkin()));

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

    for (size_t index : cellsOfInterest) {
      if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
        getCells()[index].forEach(forEachLambda, lowerCorner, higherCorner, behavior);
      }
    }
  }

  /**
   * @copydoc LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    // We increase the search region by skin, as particles can move over cell borders.
    auto startIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::subScalar(lowerCorner, this->getVerletSkin()));
    auto stopIndex3D =
        this->_cellBlock.get3DIndexOfPosition(utils::ArrayMath::addScalar(higherCorner, this->getVerletSkin()));

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

    for (size_t index : cellsOfInterest) {
      if (!_cellBlock.ignoreCellForIteration(index, behavior)) {
        getCells()[index].reduce(reduceLambda, result, lowerCorner, higherCorner, behavior);
      }
    }
  }

  /**
   * Get the cell block, not supposed to be used except by verlet lists
   * @return the cell block
   */
  internal::CellBlock3D<ReferenceCell> &getCellBlock() { return _cellBlock; }

  /**
   * @copydoc getCellBlock()
   * @note const version
   */
  const internal::CellBlock3D<ReferenceCell> &getCellBlock() const { return _cellBlock; }

  /**
   * Returns reference to the data of LinkedCellsReferences
   * @return the data
   */
  std::vector<ReferenceCell> &getCells() { return this->_cells; }

  /**
   * @copydoc getCells()
   * @note const version
   */
  const std::vector<ReferenceCell> &getCells() const { return this->_cells; }

 protected:
  /**
   * object that stores the actual Particles and keeps track of the references.
   */
  ParticleVector<ParticleType> _particleList;
  /**
   * object to manage the block of cells.
   */
  internal::CellBlock3D<ReferenceCell> _cellBlock;
  /**
   * load estimation algorithm for balanced traversals.
   */
  autopas::LoadEstimatorOption _loadEstimator;
  /**
   * Workaround for adding particles in parallel -> https://github.com/AutoPas/AutoPas/issues/555
   */
  AutoPasLock addParticleLock;
};

}  // namespace autopas
