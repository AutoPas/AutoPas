/**
 * @file HierarchicalGrid.h
 * @date 29.11.2024
 * @author atacann
 */

#pragma once

#include <algorithm>
#include <array>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/HGTraversalBase.h"
#include "autopas/containers/linkedCells/traversals/HGTraversalInterface.h"
#include "autopas/iterators/ContainerIterator.h"

namespace autopas {

/**
 * HierarchicalGrid class.
 * This class stores multiple LinkedCells containers, each for different sizes of particles
 * Traverses them all one by one, and cross hierarchy afterward
 * @tparam Particle type of the Particle
 */
template <class Particle>
class HierarchicalGrid : public ParticleContainerInterface<Particle> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle;

  /**
   * Constructor of the HierarchicalGrid class.
   * @param boxMin
   * @param boxMax
   * @param baseCutoff base cutoff for each particle, it will be scaled by the size of a Particle
   * @param cutoffs cutoffs for each level of the hierarchy
   * @param skin
   * @param rebuildFrequency
   * @param cellSizeFactor cell size factor relative to cutoff
   */
  HierarchicalGrid(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double baseCutoff,
                   const std::vector<double> &cutoffs, const double skin, const unsigned int rebuildFrequency,
                   const double cellSizeFactor = 1.0)
      : ParticleContainerInterface<Particle>(skin),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _baseCutoff(baseCutoff),
        _skin(skin),
        _numLevels(cutoffs.size()),
        _cellSizeFactor(cellSizeFactor),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(size_t)),
        _cutoffs(cutoffs) {
    /*
     * NOTE: In the future add auto cutoff? We need to know particles sizes beforehand.
     * Not relevant for md simulations but can be useful for DEM simulations.
     * Good heuristic for auto cutoffs:
     * Let m = avg particle per cell, L = number of levels
     * try to make m for each hierarchy level close, then minimize L*m
     * prob works -> ternary search for m, for fixed m binary search over particle array sorted by size to find L
     */
    // resize iteration storage variables
    _currentLevel.resize(autopas_get_max_threads() * _cacheOffset, 0);
    _prevCells.resize(autopas_get_max_threads() * _cacheOffset, 0);

    if (_cutoffs.empty()) {
      utils::ExceptionHandler::exception("Error: Hierarchical Grid cutoffs vector is empty.");
    }
    // calculate max allowable cutoff
    double maxLength = 1e15;
    for (size_t i = 0; i < 3; ++i) {
      maxLength = std::min(maxLength, boxMax[i] - boxMin[i]);
    }
    maxLength -= _skin;
    // fix cutoffs if they are bigger than maxLength
    for (auto &cutoff : _cutoffs) {
      if (cutoff > maxLength) {
        AutoPasLog(ERROR,
                   "Cutoff {} is greater than the maximum length according to Box size - skin {}, lowering cutoff",
                   cutoff, maxLength);
        cutoff = maxLength;
      }
    }
    // make sure cutoffs are sorted
    std::sort(_cutoffs.begin(), _cutoffs.end());
    // make all cutoffs unique
    _cutoffs.resize(std::unique(_cutoffs.begin(), _cutoffs.end()) - _cutoffs.begin());
    if (_cutoffs.size() != _numLevels) {
      AutoPasLog(ERROR, "Cutoffs should all be unique, adjusting cutoffs size from {} to {}", _numLevels,
                 _cutoffs.size());
      _numLevels = _cutoffs.size();
    }
    // largest interaction length
    const double interactionLengthLargest = _cutoffs.back() + _skin;
    // generate LinkedCells for each hierarchy, with different cutoffs
    _levels.reserve(_numLevels);
    for (size_t i = 0; i < _numLevels; ++i) {
      // here, LinkedCells with smaller cutoffs need to have bigger halo regions, because particles on level above
      // can potentially interact with those, even though these halo particles cannot interact with other particles
      // inside its own level a hacky way: make all LinkedCells interactionLength equal to the biggest one, adjust
      // cellSizeFactor for each LinkedCells so that the cellLength is equal to actual interaction length of that level
      const double interactionLengthLevel = _cutoffs[i] + _skin;
      const double ratio = interactionLengthLevel / interactionLengthLargest;
      _levels.emplace_back(std::make_unique<autopas::LinkedCells<Particle>>(_boxMin, _boxMax, _cutoffs.back(), skin,
                                                                            rebuildFrequency, _cellSizeFactor * ratio));
    }
  }

  /**
   * Destructor of HierarchicalGrid.
   */
  ~HierarchicalGrid() override = default;

  /**
   * Delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  HierarchicalGrid(const HierarchicalGrid &obj) = delete;

  /**
   * Delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  HierarchicalGrid &operator=(const HierarchicalGrid &other) = delete;

  /**
   * For Hierarchical Grid the cutoff is the cutoff of the level that is the largest
   * @return the largest cutoff
   */
  [[nodiscard]] double getCutoff() const final { return _cutoffs.back(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(const double cutoff) final { _baseCutoff = cutoff; }

  [[nodiscard]] double getVerletSkin() const final { return _skin; }

  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior) const override {
    size_t numParticles = 0;
    for (const auto &linkedCells : _levels) {
      numParticles += linkedCells->getNumberOfParticles(behavior);
    }
    return numParticles;
  }

  [[nodiscard]] size_t size() const override {
    size_t numParticles = 0;
    for (const auto &linkedCells : _levels) {
      numParticles += linkedCells->size();
    }
    return numParticles;
  }

  /**
   * Set the number of time-steps since last neighbor list rebuild
   * Hierarchical Grid needs to call this for each subcontainer
   * @param stepsSinceLastRebuild steps since last neighbor list rebuild
   */
  void setStepsSinceLastRebuild(size_t stepsSinceLastRebuild) override {
    AutoPasLog(INFO, "Steps since rebuild {}", stepsSinceLastRebuild);
    for (const auto &linkedCells : _levels) {
      linkedCells->setStepsSinceLastRebuild(stepsSinceLastRebuild);
    }
    ParticleContainerInterface<Particle>::setStepsSinceLastRebuild(stepsSinceLastRebuild);
  }
  /**
   * String representation of Hgrid for debugging
   * @return String showing each hierarchy and which particle ID's they contain
   */
  std::string toString() const {
    using utils::ArrayUtils::operator<<;
    std::ostringstream text;
    text << "\n------------------------------------------\nHierarchicalGrid sizes: ";
    text << utils::ArrayUtils::to_string(_cutoffs, ", ", {"Cutoffs [ ", " ]\n"});
    text << "BoxMin: " << getBoxMin() << " BoxMax: " << getBoxMax() << "\n";
    for (size_t i = 0; i < _numLevels; ++i) {
      auto &cellBlock = _levels[i]->getCellBlock();
      text << "------------------------------------------\nHierarchicalGrid " << i + 1
           << " numParticles: " << _levels[i]->size() << " cutoff: " << _cutoffs[i] << "\n";
      text << "CellLength: " << cellBlock.getCellLength()
           << " InteractionLength: " << _levels[i]->getInteractionLength()
           << "\nCellsPerInteractionLength: " << cellBlock.getCellsPerInteractionLength()
           << " NumCells: " << cellBlock.getNumCells() << "\n"
           << "CellsPerDimensionWithHalo: " << cellBlock.getCellsPerDimensionWithHalo() << "\n"
           << "Particle Density " << 1.0 * _levels[i]->size() / cellBlock.getNumCells() << "\n";
      // _levels[i]->forEach(
      //     [&](Particle &p) {
      //       text << "Cell index: " << cellBlock.get3DIndexOfPosition(p.getR()) << "\n" << p.toString() << "\n\n";
      //     },
      //     IteratorBehavior::ownedOrHaloOrDummy);
    }
    return text.str();
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const final { return _boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const final { return _boxMin; }

  /**
   * In HierarchicalGrids, interaction length returns maximum interaction length which is the biggest cutoff + skin
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const final { return _cutoffs.back() + _skin; }

  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::hierarchicalGrid; }

  [[nodiscard]] CellType getParticleCellTypeEnum() const override { return CellType::FullParticleCell; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    // reserve overestimates because reserve is designed for only one hierarchy
    // can get size() for each hierarchy and scale numParticles, numParticlesHaloEstimate proportionally?
    for (size_t i = 0; i < _numLevels; ++i) {
      _levels[i]->reserve(numParticles, numParticlesHaloEstimate);
    }
  }

  void addParticleImpl(const ParticleType &p) override { _levels[getHierarchyLevel(p)]->addParticleImpl(p); }

  bool deleteParticle(Particle &particle) override {
    return _levels[getHierarchyLevel(particle)]->deleteParticle(particle);
  }

  bool deleteParticle(size_t cellIndex, size_t particleIndex) override {
    size_t prevCells = 0;
    // find which hierarchy contains cellIndex
    for (size_t currentLevel = 0; currentLevel < _numLevels; ++currentLevel) {
      const auto cellCount = _levels[currentLevel]->getCells().size();
      if (cellIndex < prevCells + cellCount) {
        return _levels[currentLevel]->deleteParticle(cellIndex - prevCells, particleIndex);
      }
      prevCells += cellCount;
    }
    utils::ExceptionHandler::exception("Error: Couldn't find particle at cellIndex: {} particleIndex: {}", cellIndex,
                                       particleIndex);
    return false;
  }

  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    _levels[getHierarchyLevel(haloParticle)]->addHaloParticleImpl(haloParticle);
  }

  bool updateHaloParticle(const ParticleType &haloParticle) override {
    return _levels[getHierarchyLevel(haloParticle)]->updateHaloParticle(haloParticle);
  }

  void deleteHaloParticles() override {
    AUTOPAS_OPENMP(parallel for)
    for (size_t idx = 0; idx < _numLevels; idx++) {
      _levels[idx]->deleteHaloParticles();
    }
  }

  void deleteAllParticles() override {
    for (auto &linkedCells : _levels) {
      linkedCells->deleteAllParticles();
    }
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {}

  void computeInteractions(TraversalInterface *traversal) override {
    prepareTraversal(traversal);
    traversal->initTraversal();
    traversal->traverseParticles();
    traversal->endTraversal();
    // AutoPasLog(INFO, toString());
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
    std::vector<ParticleType> invalidParticles;
    // parallelized inside each updateContainer
    for (auto &linkedCells : _levels) {
      auto invalidParticlesSingle = linkedCells->updateContainer(keepNeighborListsValid);
      invalidParticles.insert(invalidParticles.end(), std::make_move_iterator(invalidParticlesSingle.begin()),
                              std::make_move_iterator(invalidParticlesSingle.end()));
    }
    return invalidParticles;
  }

  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // return traversal info of hierarchy with biggest interactionLength for autoPas container
    return _levels.back()->getTraversalSelectorInfo();
  }

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
   * @tparam regionIter
   * @param cellIndex
   * @param particleIndex
   * @param iteratorBehavior
   * @param boxMin
   * @param boxMax
   * @return tuple<ParticlePointer, CellIndex, ParticleIndex>
   */
  template <bool regionIter>
  std::tuple<const Particle *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
                                                               IteratorBehavior iteratorBehavior,
                                                               const std::array<double, 3> &boxMin,
                                                               const std::array<double, 3> &boxMax) const {
    const auto threadNum = autopas_get_thread_num();
    if (cellIndex == 0 && particleIndex == 0) {
      // new iteration
      _currentLevel[threadNum * _cacheOffset] = _prevCells[threadNum * _cacheOffset] = 0;
    }

    // get stored data
    auto &currentLevel = _currentLevel[threadNum * _cacheOffset];
    auto &prevCells = _prevCells[threadNum * _cacheOffset];

    // get which hierarchy cellIndex is in
    while (currentLevel < _numLevels) {
      const auto cellCount =
          static_cast<const CellBasedParticleContainer<FullParticleCell<Particle>> &>(*_levels[currentLevel])
              .getCells()
              .size();
      if (cellIndex - prevCells >= cellCount) {
        // skip this level
        prevCells += cellCount;
        ++currentLevel;
        continue;
      }
      const auto ret = _levels[currentLevel]->template getParticleImpl<regionIter>(cellIndex - prevCells, particleIndex,
                                                                                   iteratorBehavior, boxMin, boxMax);
      if (std::get<0>(ret) == nullptr) {
        // not found, check next level
        if (cellIndex - prevCells < cellCount) {
          // if current level was just finished, we need to increase cellIndex here
          if (iteratorBehavior & IteratorBehavior::forceSequential) {
            ++cellIndex;
          } else {
            cellIndex = prevCells + cellCount;
          }
          particleIndex = 0;
        }
        prevCells += cellCount;
        ++currentLevel;
      } else {
        return {std::get<0>(ret), std::get<1>(ret) + prevCells, std::get<2>(ret)};
      }
    }
    return {nullptr, 0, 0};
  }

  [[nodiscard]] ContainerIterator<ParticleType, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, true, true>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<ParticleType, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  [[nodiscard]] ContainerIterator<ParticleType, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<ParticleType, false, true>::ParticleVecType *additionalVectors =
          nullptr) const override {
    return ContainerIterator<ParticleType, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  [[nodiscard]] ContainerIterator<ParticleType, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, true, false>::ParticleVecType *additionalVectors = nullptr) override {
    return ContainerIterator<ParticleType, true, false>(*this, behavior, additionalVectors);
  }

  [[nodiscard]] ContainerIterator<ParticleType, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<ParticleType, false, false>::ParticleVecType *additionalVectors =
          nullptr) const override {
    return ContainerIterator<ParticleType, false, false>(*this, behavior, additionalVectors);
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
    for (size_t idx = 0; idx < _levels.size(); idx++) {
      _levels[idx]->reduce(reduceLambda, result, behavior);
    }
  }

  /**
   * Execute code on all particles in this container as defined by a lambda function.
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    for (size_t idx = 0; idx < _levels.size(); idx++) {
      _levels[idx]->forEach(forEachLambda, behavior);
    }
  }

  /**
   * Execute code on all particles in this container in a certain region as defined by a lambda function.
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    for (size_t idx = 0; idx < _levels.size(); idx++) {
      _levels[idx]->forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    }
  }

  /**
   * Execute code on all particles in this container in a certain region as defined by a lambda function.
   * @tparam Lambda (Particle &p, A &result) -> void
   * @tparam A type of reduction Value
   * @param reduceLambda code to be executed on all particles
   * @param result reference to starting and final value for reduction
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    for (size_t idx = 0; idx < _levels.size(); idx++) {
      _levels[idx]->reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    }
  }

  //
  //  /**
  // * Get the cell block, not supposed to be used except by verlet lists
  // * @return the cell block
  // */
  //  internal::CellBlock3D<ParticleCell> &getCellBlock() { return nullptr; }
  //
  //  /**
  // * @copydoc autopas::LinkedCells::getCellBlock()
  // * @note const version
  // */
  //  const internal::CellBlock3D<ParticleCell> &getCellBlock() const { return nullptr; }
  //
  //  /**
  //   * Returns a non-const reference to the cell data structure.
  //   * @return Non-const reference.
  //   */
  //  std::vector<ParticleCell> &getCells() { return nullptr; }
  //
  //  /**
  //  * Get immutable vector of cells.
  //  * @return immutable reference to _cells
  //  */
  //  [[nodiscard]] const std::vector<ParticleCell> &getCells() const {
  //    return std::vector<ParticleCell>();
  //  }

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _baseCutoff;
  double _skin;
  size_t _numLevels;
  const double _cellSizeFactor;
  // store iteration stage
  mutable std::vector<size_t> _currentLevel;
  mutable std::vector<size_t> _prevCells;
  const size_t _cacheOffset;

  /**
   *
   * @param p Particle to add into HierarchicalGrid
   * @return which Hierarchy the particle belongs to
   */
  size_t getHierarchyLevel(const ParticleType &p) const {
    // scale size by baseCutoff
    const double cutoff = p.getSize() * _baseCutoff;
    // AutoPasLog(INFO, "Adding particle with size {}\n{}", ParticleType::getSizeFunctor(&p), p.toString());
    for (size_t i = 0; i < _numLevels; ++i) {
      if (_cutoffs[i] >= cutoff) {
        return i;
      }
    }
    AutoPasLog(ERROR,
               "Size of Particle times baseCutoff is bigger than biggest cutoff of HierarchicalGrid, "
               "will result in wrong interaction calculation.");
    return _numLevels - 1;
  }

 protected:
  /**
   * Checks if a given traversal is allowed for HierarchicalGrid and sends required data to HGridTraversal
   * @param traversal
   */
  void prepareTraversal(TraversalInterface *traversal) {
    auto *traversalInterface = dynamic_cast<HGTraversalInterface *>(traversal);
    auto *hGridTraversal = dynamic_cast<HGTraversalBase<ParticleCell> *>(traversal);
    if (traversalInterface && hGridTraversal) {
      hGridTraversal->setLevels(&_levels, _cutoffs, _skin, this->getMaxDisplacement());
    } else {
      autopas::utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the HierarchicalGrid container. TraversalID: {}",
          traversal->getTraversalType());
    }
  }

  std::vector<std::unique_ptr<LinkedCells<Particle>>> _levels;
  std::vector<double> _cutoffs;
};

}  // namespace autopas