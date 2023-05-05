/**
 * @file Octree.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#pragma once

#include <cstdio>
#include <list>
#include <stack>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/containers/LeavingParticleCollector.h"
#include "autopas/containers/octree/OctreeLeafNode.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeNodeWrapper.h"
#include "autopas/containers/octree/traversals/OTTraversalInterface.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/utils/ParticleCellHelpers.h"
#include "autopas/utils/inBox.h"
#include "autopas/utils/logging/OctreeLogger.h"

namespace autopas {

/**
 * The octree is a CellBasedParticleContainer that consists internally of two octrees. One for owned and one for halo
 * particles. It abuses the CellBasedParticleContainer::_cell vector to hold the root nodes for each.
 *
 * The tree consists of OctreeNodeWrapper objects, which
 *
 * @note Octree has a particular to interpret the index of the ContainerIterator. For details see getParticleImpl().
 *
 * @tparam Particle
 */
template <class Particle>
class Octree : public CellBasedParticleContainer<OctreeNodeWrapper<Particle>>,
               public internal::CellBorderAndFlagManager {
 public:
  /**
   * The particle cell used in this CellBasedParticleContainer
   */
  using ParticleCell = OctreeNodeWrapper<Particle>;

  /**
   * The particle type used in this container.
   */
  using ParticleType = typename ParticleCell::ParticleType;

  /**
   * This particle container contains two cells. Both cells are octrees.
   * - this->_cells[CellTypes::OWNED] yields the octree that contains all owned particles
   * - this->_cells[CellTypes::HALO] yields the octree that contains all halo particles
   */
  enum CellTypes : int { OWNED = 0, HALO = 1 };

  /**
   * A cell index that is definitely always invalid.
   */
  constexpr static size_t invalidCellIndex = 9;

  /**
   * Construct a new octree with two sub-octrees: One for the owned particles and one for the halo particles.
   * @param boxMin The minimum coordinate of the enclosing box
   * @param boxMax The maximum coordinate of the enclosing box
   * @param cutoff The cutoff radius
   * @param skinPerTimestep The skin radius per timestep
   * @param rebuildFrequency The Rebuild Frequency
   * @param cellSizeFactor The cell size factor
   */
  Octree(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double cutoff,
         const double skinPerTimestep, const unsigned int rebuildFrequency, const double cellSizeFactor)
      : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skinPerTimestep * rebuildFrequency) {
    using namespace autopas::utils::ArrayMath::literals;

    // @todo Obtain this from a configuration, reported in https://github.com/AutoPas/AutoPas/issues/624
    int unsigned treeSplitThreshold = 16;

    double interactionLength = this->getInteractionLength();

    // Create the octree for the owned particles
    this->_cells.push_back(
        OctreeNodeWrapper<Particle>(boxMin, boxMax, treeSplitThreshold, interactionLength, cellSizeFactor));

    // Extend the halo region with cutoff + skin in all dimensions
    auto haloBoxMin = boxMin - interactionLength;
    auto haloBoxMax = boxMax + interactionLength;
    // Create the octree for the halo particles
    this->_cells.push_back(
        OctreeNodeWrapper<Particle>(haloBoxMin, haloBoxMax, treeSplitThreshold, interactionLength, cellSizeFactor));
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer(bool keepNeighborListValid) override {
    // invalidParticles: all outside boxMin/Max
    std::vector<Particle> invalidParticles{};

    if (keepNeighborListValid) {
      invalidParticles = LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    } else {
      // This is a very primitive and inefficient way to rebuild the container:

      // @todo Make this less indirect. (Find a better way to iterate all particles inside the octree to change
      //   this function back to a function that actually copies all particles out of the octree.)
      //   The problem is captured by https://github.com/AutoPas/AutoPas/issues/622

      // 1. Copy all particles out of the container
      std::vector<Particle *> particleRefs;
      this->_cells[CellTypes::OWNED].collectAllParticles(particleRefs);
      std::vector<Particle> particles{};
      particles.reserve(particleRefs.size());

      for (auto *p : particleRefs) {
        if (p->isDummy()) {
          // don't do anything with dummies. They will just be dropped when the container is rebuilt.
          continue;
        } else if (utils::inBox(p->getR(), this->getBoxMin(), this->getBoxMax())) {
          particles.push_back(*p);
        } else {
          invalidParticles.push_back(*p);
        }
      }

      // 2. Clear the container
      this->deleteAllParticles();

      // 3. Insert the particles back into the container
      for (auto &particle : particles) {
        addParticleImpl(particle);
      }
    }

    return invalidParticles;
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    if (auto *traversalInterface = dynamic_cast<OTTraversalInterface<ParticleCell> *>(traversal)) {
      traversalInterface->setCells(&this->_cells);
    }

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::octree; }

  /**
   * @copydoc ParticleContainerInterface::getParticleCellTypeEnum()
   * @note The octree stores the particles in its leaves. Even though the leaves inherit from `FullParticleCell`, they
   * are an extension of `FullParticleCell`. However, the interface stays the same and the leaves can be treated just
   * like regular `FullParticleCell`s.
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() override { return CellType::FullParticleCell; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    // TODO create a balanced tree and reserve space in the leaves.
  }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const ParticleType &p) override { this->_cells[CellTypes::OWNED].addParticle(p); }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    ParticleType p_copy = haloParticle;
    p_copy.setOwnershipState(OwnershipState::halo);
    this->_cells[CellTypes::HALO].addParticle(p_copy);
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    ParticleType pCopy = haloParticle;
    pCopy.setOwnershipState(OwnershipState::halo);
    return internal::checkParticleInCellAndUpdateByIDAndPosition(this->_cells[CellTypes::HALO], pCopy,
                                                                 this->getVerletSkin());
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {}

  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                           IteratorBehavior iteratorBehavior,
                                                           const std::array<double, 3> &boxMin,
                                                           const std::array<double, 3> &boxMax) const override {
    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }
  std::tuple<const Particle *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                           IteratorBehavior iteratorBehavior) const override {
    // this is not a region iter hence we stretch the bounding box to the numeric max
    constexpr std::array<double, 3> boxMin{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(),
                                           std::numeric_limits<double>::lowest()};

    constexpr std::array<double, 3> boxMax{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                           std::numeric_limits<double>::max()};
    return getParticleImpl<false>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }

  /**
   * Container specific implementation for getParticle. See ParticleContainerInterface::getParticle().
   *
   * @note In this context cell == leaf cell
   * @note The index encodes the location in the octree. Each digit signifies which child of the node to enter.
   * The right most digit selects the tree, the next the child of the root, and so on. So for example 31 would be the
   * fourth child (3) in the halo tree (1). If this is not a leaf node, we recursively follow the first children
   * (= prepend zeros) until the deepest level is found.
   *
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @return tuple<ParticlePointer, CellIndex, ParticleIndex>
   */
  template <bool regionIter>
  std::tuple<const ParticleType *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
                                                                   IteratorBehavior iteratorBehavior,
                                                                   const std::array<double, 3> &boxMin,
                                                                   const std::array<double, 3> &boxMax) const {
    // FIXME think about parallelism.
    // This `if` currently disables it but should be replaced with logic that determines the start index.
    if (autopas_get_thread_num() > 0 and not(iteratorBehavior & IteratorBehavior::forceSequential)) {
      return {nullptr, 0, 0};
    }
    // if owned particles are not interesting jump directly to the halo tree
    // FIXME: with num threads > 1 the start cell IDs will have to be set to more complicated values here
    if (cellIndex == 0 and not(iteratorBehavior & IteratorBehavior::owned)) {
      cellIndex = HALO;
    }

    // shortcut if the given index doesn't exist
    if (cellIndex < 10 and cellIndex > HALO) {
      return {nullptr, 0, 0};
    }

    std::vector<size_t> currentCellIndex{};
    OctreeLeafNode<Particle> *currentCellPtr = nullptr;

    std::tie(currentCellIndex, currentCellPtr) = getLeafCellByIndex(cellIndex);
    // check the data behind the indices
    if (particleIndex >= currentCellPtr->getNumberOfParticles() or
        not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
            (*currentCellPtr)[particleIndex], iteratorBehavior, boxMin, boxMax)) {
      // either advance them to something interesting or invalidate them.
      std::tie(currentCellPtr, particleIndex) = advanceIteratorIndices<regionIter>(
          currentCellIndex, currentCellPtr, particleIndex, iteratorBehavior, boxMin, boxMax);
    }

    // shortcut if the given index doesn't exist
    if (currentCellPtr == nullptr) {
      return {nullptr, 0, 0};
    }
    // parse cellIndex and get referenced cell and particle
    const Particle *retPtr = &((*currentCellPtr)[particleIndex]);

    // if no err value was set, convert cell index from vec to integer
    if (currentCellIndex.empty()) {
      cellIndex = invalidCellIndex;
    } else {
      cellIndex = 0;
      // needs signed int to prevent underflow
      for (int i = static_cast<int>(currentCellIndex.size()) - 1; i >= 0; --i) {
        cellIndex *= 10;
        cellIndex += currentCellIndex[i];
      }
    }

    return {retPtr, cellIndex, particleIndex};
  }

  /**
   * Helper function to retrieve the pointer to a leaf cell as well as properly parse the cell index.
   *
   * @note This function assumes the given index is valid. If it is not some exceptions might be thrown or operator[]
   * might access invalid memory
   *
   * @param cellIndex Index in the tree as integer
   * @return tuple<index as vector, pointer to cell>
   */
  std::tuple<std::vector<size_t>, OctreeLeafNode<Particle> *> getLeafCellByIndex(size_t cellIndex) const {
    // parse cellIndex and get referenced cell and particle
    std::vector<size_t> currentCellIndex;
    // constant heuristic for the tree depth
    currentCellIndex.reserve(10);
    currentCellIndex.push_back(cellIndex % 10);
    cellIndex /= 10;
    OctreeNodeInterface<Particle> *currentCell = this->_cells[currentCellIndex.back()].getRaw();
    // don't restrict loop via cellIndex because it might have "hidden" leading 0
    while (currentCell->hasChildren()) {
      currentCellIndex.push_back(cellIndex % 10);
      cellIndex /= 10;
      currentCell = currentCell->getChild(currentCellIndex.back());
    }
    return {currentCellIndex, dynamic_cast<OctreeLeafNode<Particle> *>(currentCell)};
  }

  bool deleteParticle(Particle &particle) override {
    if (particle.isOwned()) {
      return this->_cells[CellTypes::OWNED].deleteParticle(particle);
    } else if (particle.isHalo()) {
      return this->_cells[CellTypes::HALO].deleteParticle(particle);
    } else {
      utils::ExceptionHandler::exception("Particle to be deleted is neither owned nor halo!\n" + particle.toString());
      return false;
    }
  }

  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {
    auto [cellIndexVector, cell] = getLeafCellByIndex(cellIndex);
    auto &particleVec = cell->_particles;
    auto &particle = particleVec[particleIndex];
    // swap-delete
    particle = particleVec.back();
    particleVec.pop_back();
    return particleIndex < particleVec.size();
  }

  [[nodiscard]] ContainerIterator<Particle, true, false> begin(
      IteratorBehavior behavior,
      typename ContainerIterator<Particle, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<Particle, true, false>(*this, behavior, additionalVectors);
  }

  [[nodiscard]] ContainerIterator<Particle, false, false> begin(
      IteratorBehavior behavior,
      typename ContainerIterator<Particle, false, false>::ParticleVecType *additionalVectors = nullptr) const override {
    return ContainerIterator<Particle, false, false>(*this, behavior, additionalVectors);
  }

  [[nodiscard]] ContainerIterator<Particle, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle, true, true>::ParticleVecType *additionalVectors) override {
    return ContainerIterator<Particle, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  [[nodiscard]] ContainerIterator<Particle, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle, false, true>::ParticleVecType *additionalVectors) const override {
    return ContainerIterator<Particle, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    using namespace autopas::utils::ArrayMath::literals;

    // this is a dummy since it is not actually used
    std::array<unsigned long, 3> dims = {1, 1, 1};
    std::array<double, 3> cellLength = this->getBoxMax() - this->getBoxMin();
    return TraversalSelectorInfo(dims, this->getInteractionLength(), cellLength, 0);
  }

  /**
   * Get the number of particles that belong to this octree. (Owned and and halo.)
   *
   * @return The integer # of particles in the container
   */
  [[nodiscard]] unsigned long getNumberOfParticles() const override {
    return this->_cells[CellTypes::OWNED].numParticles() + this->_cells[CellTypes::HALO].numParticles();
  }

  void deleteHaloParticles() override { this->_cells[CellTypes::HALO].clear(); }

  [[nodiscard]] bool cellCanContainHaloParticles(std::size_t i) const override {
    if (i > 1) {
      throw std::runtime_error("[Octree.h]: This cell container (octree) contains only two cells");
    }
    return i == CellTypes::HALO;
  }

  [[nodiscard]] bool cellCanContainOwnedParticles(std::size_t i) const override {
    if (i > 1) {
      throw std::runtime_error("[Octree.h]: This cell container (octree) contains only two cells");
    }
    return i == CellTypes::OWNED;
  }

  /**
   * Execute code on all particles in this container as defined by a lambda function.
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior & IteratorBehavior::owned) this->_cells[OWNED].forEach(forEachLambda);
    if (behavior & IteratorBehavior::halo) this->_cells[HALO].forEach(forEachLambda);
    if (not(behavior & IteratorBehavior::ownedOrHalo))
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownedOrHalo
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    if (behavior & IteratorBehavior::owned) this->_cells[OWNED].reduce(reduceLambda, result);
    if (behavior & IteratorBehavior::halo) this->_cells[HALO].reduce(reduceLambda, result);
    if (not(behavior & IteratorBehavior::ownedOrHalo))
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
  }

  /**
   * @copydoc LinkedCells::forEachInRegion()
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned)
      this->_cells[OWNED].forEachInRegion(forEachLambda, lowerCorner, higherCorner);
    if (behavior & IteratorBehavior::halo) this->_cells[HALO].forEachInRegion(forEachLambda, lowerCorner, higherCorner);
    if (not(behavior & IteratorBehavior::ownedOrHalo))
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
  }

  /**
   * @copydoc LinkedCells::reduceInRegion()
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    if (behavior & IteratorBehavior::owned)
      this->_cells[OWNED].reduceInRegion(reduceLambda, result, lowerCorner, higherCorner);
    if (behavior & IteratorBehavior::halo)
      this->_cells[HALO].reduceInRegion(reduceLambda, result, lowerCorner, higherCorner);
    if (not(behavior & IteratorBehavior::ownedOrHalo))
      utils::ExceptionHandler::exception("Encountered invalid iterator behavior!");
  }

 private:
  /**
   * Given a pair of cell-/particleIndex and iterator restrictions either returns the next indices that match these
   * restrictions or indices that are out of bounds.
   *
   * @tparam regionIter
   * @param currentCellIndex In/Out parameter
   * @param currentCellPtr
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @return tuple<nextLeafCell, nextParticleIndex>
   */
  template <bool regionIter>
  std::tuple<OctreeLeafNode<Particle> *, size_t> advanceIteratorIndices(
      std::vector<size_t> &currentCellIndex, OctreeNodeInterface<Particle> *const currentCellPtr, size_t particleIndex,
      IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin,
      const std::array<double, 3> &boxMax) const {
    // TODO: parallelize at the higher tree levels. Choose tree level to parallelize via log_8(numThreads)
    const size_t minLevel = 0;
    //        (iteratorBehavior & IteratorBehavior::forceSequential) or autopas_get_num_threads() == 1
    //            ? 0
    //            : static_cast<size_t>(std::ceil(std::log(static_cast<double>(autopas_get_num_threads())) /
    //            std::log(8.)));
    OctreeNodeInterface<Particle> *currentCellInterfacePtr = currentCellPtr;
    OctreeLeafNode<Particle> *currentLeafCellPtr = nullptr;

    // helper function:
    auto cellIsRelevant = [&](const OctreeNodeInterface<Particle> *const cellPtr) {
      bool isRelevant = cellPtr->getNumberOfParticles() > 0;
      if constexpr (regionIter) {
        // particles can move over cell borders. Calculate the volume this cell's particles can be.
        const auto cellLowCornerSkin = utils::ArrayMath::subScalar(cellPtr->getBoxMin(), this->getVerletSkin() * 0.5);
        const auto cellHighCornerSkin = utils::ArrayMath::addScalar(cellPtr->getBoxMax(), this->getVerletSkin() * 0.5);
        isRelevant = utils::boxesOverlap(cellLowCornerSkin, cellHighCornerSkin, boxMin, boxMax);
      }
      return isRelevant;
    };

    // find the next particle of interest which might be in a different cell
    do {
      // advance to the next particle
      ++particleIndex;

      // this loop finds the next relevant leaf cell or triggers a return
      // flag for weird corner cases. See further down.
      bool forceJumpToNextCell = false;
      // If this breaches the end of a cell, find the next non-empty cell and reset particleIndex.
      while (particleIndex >= currentCellInterfacePtr->getNumberOfParticles() or forceJumpToNextCell) {
        // CASE: we are at the end of a branch
        // => Move up until there are siblings that were not touched yet
        while (currentCellIndex.back() == 7) {
          currentCellInterfacePtr = currentCellInterfacePtr->getParent();
          currentCellIndex.pop_back();
          // If there are no more cells in the branch that we are responsible for set invalid parameters and return.
          if (currentCellIndex.size() < minLevel or currentCellIndex.empty()) {
            currentCellIndex.clear();
            return {nullptr, std::numeric_limits<decltype(particleIndex)>::max()};
          }
        }

        // Switch to the next child of the same parent
        ++currentCellIndex.back();
        // Identify the next (inner) cell pointer
        if (currentCellIndex.size() == 1) {
          // if we are already beyond everything
          if (currentCellIndex.back() > HALO) {
            return {nullptr, std::numeric_limits<decltype(particleIndex)>::max()};
          }
          // special case: the current thread should ALSO iterate halo particles
          if ((iteratorBehavior & IteratorBehavior::halo)
              /* FIXME: for parallelization: and ((iteratorBehavior & IteratorBehavior::forceSequential) or autopas_get_num_threads() == 1) */) {
            currentCellInterfacePtr = this->_cells[HALO].getRaw();
          } else {
            // don't jump to the halo tree -> set invalid parameters and return
            currentCellIndex.clear();
            return {nullptr, std::numeric_limits<decltype(particleIndex)>::max()};
          }
        } else {
          currentCellInterfacePtr = currentCellInterfacePtr->getParent()->getChild(currentCellIndex.back());
        }
        // check that the inner cell is actually interesting, otherwise skip it.
        if (not cellIsRelevant(currentCellInterfacePtr)) {
          forceJumpToNextCell = true;
          continue;
        }
        // The inner cell is relevant so descend to its first relevant leaf
        forceJumpToNextCell = false;
        while (currentCellInterfacePtr->hasChildren() and not forceJumpToNextCell) {
          // find the first child in the relevant region
          size_t firstRelevantChild = 0;
          // if the child is irrelevant (empty or outside the region) skip it
          const auto *childToCheck = currentCellInterfacePtr->getChild(firstRelevantChild);
          while (not cellIsRelevant(childToCheck)) {
            ++firstRelevantChild;
            childToCheck = currentCellInterfacePtr->getChild(firstRelevantChild);
            if (firstRelevantChild > 7) {
              // weird corner case: we descended into this branch because it overlaps with the region and has particles
              // BUT the particles are not inside the children which are in the region,
              // hence no child fulfills both requirements and all are irrelevant.
              forceJumpToNextCell = true;
              break;
            }
          }
          currentCellIndex.push_back(firstRelevantChild);
          currentCellInterfacePtr = currentCellInterfacePtr->getChild(firstRelevantChild);
        }
        particleIndex = 0;
      }

      // at this point we should point to a leaf. All other cases should have hit a return earlier.
      currentLeafCellPtr = dynamic_cast<OctreeLeafNode<Particle> *>(currentCellInterfacePtr);
      // sanity check
      if (currentLeafCellPtr == nullptr) {
        utils::ExceptionHandler::exception("Expected a leaf node but didn't get one!");
      }

    } while (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
        (*currentLeafCellPtr)[particleIndex], iteratorBehavior, boxMin, boxMax));
    return {currentLeafCellPtr, particleIndex};
  }

  /**
   * A logger that can be called to log the octree data structure.
   */
  OctreeLogger<Particle> logger;

  double skin;
};

}  // namespace autopas
