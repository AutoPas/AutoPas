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
#include "autopas/containers/hierarchicalGrid/traversals/HGTraversalBase.h"
#include "autopas/containers/hierarchicalGrid/traversals/HGTraversalInterface.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/iterators/ContainerIterator.h"

namespace autopas {

/**
 * HierarchicalGrid class.
 * This class stores multiple LinkedCells containers, each for different sizes of particles
 * Traverses them all one by one, and cross hierarchy afterward
 * @tparam Particle_T type of the Particle
 */
template <class Particle_T>
class HierarchicalGrid : public ParticleContainerInterface<Particle_T> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle_T>;
  /**
   * Constructor of the HierarchicalGrid class.
   * @param boxMin
   * @param boxMax
   * @param hGridMaxCutoffPerLevel Max cutoffs for each level of the hierarchy
   * @param skin Verlet skin
   * @param cellSizeFactor Cell size factor relative to cutoff
   * @param sortingThreshold Number of particles in two cells from which sorting should be performed
   * @param loadEstimator The load estimation algorithm for balanced traversals.
   * @param fittedGrids Whether to fit the grids into one another.
   */
  HierarchicalGrid(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                   const std::vector<double> &hGridMaxCutoffPerLevel, const double skin,
                   const double cellSizeFactor = 1.0, const size_t sortingThreshold = 8,
                   LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell,
                   bool fittedGrids = false)
      : ParticleContainerInterface<Particle_T>(skin),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _skin(skin),
        _numLevels(hGridMaxCutoffPerLevel.size()),
        _cellSizeFactor(cellSizeFactor),
        _cacheOffset(DEFAULT_CACHE_LINE_SIZE / sizeof(size_t)),
        _maxCutoffPerLevel(hGridMaxCutoffPerLevel),
        _fittedGrids(fittedGrids) {
    /*
     * @todo: In the future automatically choose the number of levels and min cell sizes per level?
     * We need to know particles sizes beforehand.
     * Not relevant for md simulations but can be useful for DEM simulations.
     * Good heuristic for automatically choosing:
     * Let m = avg particle per cell, L = number of levels
     * try to make m for each hierarchy level close, then minimize L*m
     * prob works -> ternary search for m, for fixed m binary search over particle array sorted by size to find L
     */
    // resize iteration storage variables
    _currentLevel.resize(autopas_get_max_threads() * _cacheOffset, 0);
    _prevCells.resize(autopas_get_max_threads() * _cacheOffset, 0);

    // In case cutoffs are empty, the traversals will return not applicable, as it makes no sense to use HG, instead of
    // LC

    // calculate max allowable cutoff
    double maxLength = 1e15;
    for (size_t i = 0; i < 3; ++i) {
      maxLength = std::min(maxLength, boxMax[i] - boxMin[i]);
    }
    maxLength -= _skin;
    // exception if max Cutoffs are bigger than maxLength
    for (auto &cellSize : _maxCutoffPerLevel) {
      if (cellSize > maxLength) {
        utils::ExceptionHandler::exception(
            "Cell size {} is greater than the maximum length according to Box size - skin {}", cellSize, maxLength);
      }
    }
    // make sure maxCutoffs are sorted
    std::sort(_maxCutoffPerLevel.begin(), _maxCutoffPerLevel.end());
    // make all maxCutoffs unique
    _maxCutoffPerLevel.resize(std::unique(_maxCutoffPerLevel.begin(), _maxCutoffPerLevel.end()) -
                              _maxCutoffPerLevel.begin());
    if (_maxCutoffPerLevel.size() != _numLevels) {
      AutoPasLog(WARN, "Max cutoffs should all be unique, adjusting max cutoffs size from {} to {}", _numLevels,
                 _maxCutoffPerLevel.size());
      _numLevels = _maxCutoffPerLevel.size();
    }
    if (_fittedGrids) {
      createFittedGrids(sortingThreshold, loadEstimator);
      return;
    }
    // largest interaction length
    const double interactionLengthLargest = _maxCutoffPerLevel.back() + _skin;
    // generate LinkedCells for each hierarchy, with different cutoffs
    _levels.reserve(_numLevels);
    for (size_t i = 0; i < _numLevels; ++i) {
      // here, LinkedCells with smaller cellSizes need to have bigger halo regions, because particles on level above
      // can potentially interact with those, even though these halo particles cannot interact with other particles
      // inside its own level a hacky way: make all LinkedCells interactionLength equal to the biggest one, adjust
      // cellSizeFactor for each LinkedCells so that the cellLength is equal to actual interaction length of that level
      // For traversals, we also choose the proper, smaller interaction length
      const double interactionLengthLevel = _maxCutoffPerLevel[i] + _skin;
      const double ratio = interactionLengthLevel / interactionLengthLargest;
      _levels.emplace_back(
          std::make_unique<autopas::LinkedCells<Particle_T>>(_boxMin, _boxMax, _maxCutoffPerLevel.back(), _skin,
                                                             _cellSizeFactor * ratio, sortingThreshold, loadEstimator));
    }
  }

 private:
  /**
   * Create hierarchy levels whose cell grids are fitted to each other.
   *
   * @param sortingThreshold Number of particles in two cells from which sorting should be performed.
   * @param loadEstimator The load estimation algorithm for balanced traversals.
   */
  void createFittedGrids(const size_t sortingThreshold, LoadEstimatorOption loadEstimator) {
    std::vector<double> maxInterLenPerLevel(_numLevels);
    for (size_t i = 0; i < _numLevels; i++) {
      maxInterLenPerLevel[i] = _maxCutoffPerLevel[i] + _skin;
    }

    // Assure at least 2 levels by increasing the cell length of the largest level if necessary, to make sure the
    // smallest level can fit at least 2 cells into one cell of the largest level.
    double highLevelRatio = 1.0;
    if (2.0 * maxInterLenPerLevel[0] >= maxInterLenPerLevel.back()) {
      // Smallest representable value above the exact threshold that still guarantees >=2 lower cells per upper cell.
      const double minimumSafeRatio = 2.0 * maxInterLenPerLevel[0] / maxInterLenPerLevel.back();
      highLevelRatio = std::nextafter(minimumSafeRatio, std::numeric_limits<double>::infinity());
      AutoPasLog(DEBUG,
                 "Adjusted high-level cell size by factor {:.6g} to guarantee at least two lower cells per "
                 "upper cell (smallest interaction length: {:.6g}, largest interaction length: {:.6g}).",
                 highLevelRatio, maxInterLenPerLevel[0], maxInterLenPerLevel.back());
    }

    _levels.resize(_numLevels);

    // construct last/largest level first, so we can fit other levels to it
    _levels.at(_numLevels - 1) = std::make_unique<autopas::LinkedCells<Particle_T>>(
        _boxMin, _boxMax, _maxCutoffPerLevel.back(), _skin, _cellSizeFactor * highLevelRatio, sortingThreshold,
        loadEstimator);

    auto highestCellLength = _levels.at(_numLevels - 1)->getCellBlock().getCellLength();
    std::array<size_t, 3> cellsPerDimensionHigher = {
        static_cast<size_t>(std::llround((_boxMax[0] - _boxMin[0]) / highestCellLength[0])),
        static_cast<size_t>(std::llround((_boxMax[1] - _boxMin[1]) / highestCellLength[1])),
        static_cast<size_t>(std::llround((_boxMax[2] - _boxMin[2]) / highestCellLength[2]))};
    std::array<size_t, 3> cellsPerDimension;
    size_t numLowerCellsPerHigher;

    // generate LinkedCells for each level with different cutoffs
    for (size_t i = _numLevels - 1; i-- > 0;) {
      // increase size of lower level cells to fit the bigger level cell, if it does not fit, remove that level
      auto cellLength = _levels[i + 1]->getCellBlock().getCellLength();
      for (size_t d = 0; d < 3; ++d) {
        numLowerCellsPerHigher =
            static_cast<size_t>(std::floor(cellLength[d] / (maxInterLenPerLevel[i] * _cellSizeFactor)));
        // sort out levels too close to be fitted
        if (numLowerCellsPerHigher < 2) {
          _maxCutoffPerLevel.erase(_maxCutoffPerLevel.begin() + i);
          maxInterLenPerLevel.erase(maxInterLenPerLevel.begin() + i);
          _levels.erase(_levels.begin() + i);
          _numLevels--;
          break;
        }
        // update cellsPerDimension for the next lower level
        cellsPerDimension[d] = numLowerCellsPerHigher * cellsPerDimensionHigher[d];
      }
      // if level was removed, create no linked cells and continue with next level
      if (numLowerCellsPerHigher < 2) {
        continue;
      }
      cellsPerDimensionHigher = cellsPerDimension;
      // @todo: This case is currently dead code and can be removed, I left it for now, so there is some trace of this
      // existing, in case we might want non-fitted Hgrids using the HGC08 traversal in the future.
      // Case we dont want to fit the grids perfecty, but just want to make sure they are not too close together, to
      // keep stride of 3
      if (!_fittedGrids) {
        double ratio = maxInterLenPerLevel[i] / maxInterLenPerLevel.back();
        // again the hacky way to get the proper halo size
        _levels[i] = std::make_unique<autopas::LinkedCells<Particle_T>>(_boxMin, _boxMax, _maxCutoffPerLevel.back(),
                                                                        _skin, _cellSizeFactor * ratio,
                                                                        sortingThreshold, loadEstimator);
      } else {
        // again the hacky way to get the proper halo size
        _levels[i] = std::make_unique<autopas::LinkedCells<Particle_T>>(
            _boxMin, _boxMax, _maxCutoffPerLevel.back(), _skin, cellsPerDimension, sortingThreshold, loadEstimator);
      }
    }
    _levels.shrink_to_fit();
    _maxCutoffPerLevel.shrink_to_fit();
  }

 public:
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
  [[nodiscard]] double getCutoff() const final { return _maxCutoffPerLevel.back(); }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(const double cutoff) final {}

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
    for (const auto &linkedCells : _levels) {
      linkedCells->setStepsSinceLastRebuild(stepsSinceLastRebuild);
    }
    ParticleContainerInterface<Particle_T>::setStepsSinceLastRebuild(stepsSinceLastRebuild);
  }
  /**
   * String representation of Hgrid for debugging
   * @return String showing each hierarchy and which particle ID's they contain
   */
  std::string toString() const {
    using utils::ArrayUtils::operator<<;
    std::ostringstream text;
    text << "\n------------------------------------------\nHierarchicalGrid sizes: ";
    text << utils::ArrayUtils::to_string(_maxCutoffPerLevel, ", ", {"Cutoffs [ ", " ]\n"});
    text << "BoxMin: " << getBoxMin() << " BoxMax: " << getBoxMax() << "\n";
    for (size_t i = 0; i < _numLevels; ++i) {
      auto &cellBlock = _levels[i]->getCellBlock();
      text << "------------------------------------------\nHierarchicalGrid " << i + 1
           << " numParticles: " << _levels[i]->size() << " cutoff: " << _maxCutoffPerLevel[i] << "\n";
      text << "CellLength: " << cellBlock.getCellLength()
           << " InteractionLength: " << _levels[i]->getInteractionLength()
           << "\nCellsPerInteractionLength: " << cellBlock.getCellsPerInteractionLength()
           << " NumCells: " << cellBlock.getNumCells() << "\n"
           << "CellsPerDimensionWithHalo: " << cellBlock.getCellsPerDimensionWithHalo() << "\n"
           << "Particle Density " << 1.0 * _levels[i]->size() / cellBlock.getNumCells() << "\n";
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
  [[nodiscard]] double getInteractionLength() const final { return _maxCutoffPerLevel.back() + _skin; }

  [[nodiscard]] ContainerOption getContainerType() const override {
    if (_fittedGrids) {
      return ContainerOption::hierarchicalGridFitted;
    } else {
      return ContainerOption::hierarchicalGrid;
    }
  }

  /**
   * Get the cell type used by this container.
   * @return CellType::FullParticleCell.
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() const { return CellType::FullParticleCell; }

  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    // reserve on each level proportionally to their number of cells
    double total = 0;
    for (size_t i = 0; i < _numLevels; ++i) {
      const auto cellsPerDim = _levels[i]->getTraversalSelectorInfo().cellsPerDim;
      total += 1.0 * cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2];
    }
    for (size_t i = 0; i < _numLevels; ++i) {
      const auto cellsPerDim = _levels[i]->getTraversalSelectorInfo().cellsPerDim;
      size_t numParticlesEstimatePerLevel =
          numParticles * (1.0 * cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2]) / total;
      size_t numHaloParticlesEstimatePerLevel =
          numParticlesHaloEstimate * (1.0 * cellsPerDim[0] * cellsPerDim[1] * cellsPerDim[2]) / total;
      _levels[i]->reserve(numParticlesEstimatePerLevel, numHaloParticlesEstimatePerLevel);
    }
  }

  void addParticleImpl(const Particle_T &p) override { _levels[getHierarchyLevel(p)]->addParticleImpl(p); }

  bool deleteParticle(Particle_T &particle) override {
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

  void addHaloParticleImpl(const Particle_T &haloParticle) override {
    _levels[getHierarchyLevel(haloParticle)]->addHaloParticleImpl(haloParticle);
  }

  bool updateHaloParticle(const Particle_T &haloParticle) override {
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

  [[nodiscard]] std::vector<Particle_T> updateContainer(bool keepNeighborListsValid) override {
    if (keepNeighborListsValid) {
      return autopas::LeavingParticleCollector::collectParticlesAndMarkNonOwnedAsDummy(*this);
    }
    std::vector<Particle_T> invalidParticles;
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
    auto infoHighest = _levels.back()->getTraversalSelectorInfo();
    return TraversalSelectorInfo(infoHighest.cellsPerDim, infoHighest.interactionLength, infoHighest.cellLength, 0,
                                 _numLevels);
  }

  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                             IteratorBehavior iteratorBehavior,
                                                             const std::array<double, 3> &boxMin,
                                                             const std::array<double, 3> &boxMax) const override {
    return getParticleImpl<true>(cellIndex, particleIndex, iteratorBehavior, boxMin, boxMax);
  }

  std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
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
  std::tuple<const Particle_T *, size_t, size_t> getParticleImpl(size_t cellIndex, size_t particleIndex,
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
      const auto cellCount = _levels[currentLevel]->getCells().size();
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

  [[nodiscard]] ContainerIterator<Particle_T, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      utils::optRef<typename ContainerIterator<Particle_T, true, true>::ParticleVecType> additionalVectors =
          std::nullopt) override {
    return ContainerIterator<Particle_T, true, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  [[nodiscard]] ContainerIterator<Particle_T, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      utils::optRef<typename ContainerIterator<Particle_T, false, true>::ParticleVecType> additionalVectors =
          std::nullopt) const override {
    return ContainerIterator<Particle_T, false, true>(*this, behavior, additionalVectors, lowerCorner, higherCorner);
  }

  [[nodiscard]] ContainerIterator<Particle_T, true, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      utils::optRef<typename ContainerIterator<Particle_T, true, false>::ParticleVecType> additionalVectors =
          std::nullopt) override {
    return ContainerIterator<Particle_T, true, false>(*this, behavior, additionalVectors);
  }

  [[nodiscard]] ContainerIterator<Particle_T, false, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo,
      utils::optRef<typename ContainerIterator<Particle_T, false, false>::ParticleVecType> additionalVectors =
          std::nullopt) const override {
    return ContainerIterator<Particle_T, false, false>(*this, behavior, additionalVectors);
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle_T p, A initialValue) -> void
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
   * @tparam Lambda (Particle_T &p) -> void
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
   * @tparam Lambda (Particle_T &p) -> void
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
   * @tparam Lambda (Particle_T &p, A &result) -> void
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

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _skin;
  size_t _numLevels;
  const double _cellSizeFactor;
  /**
   * Per-thread index of the hierarchy level currently used by getParticleImpl().
   */
  mutable std::vector<size_t> _currentLevel;
  /**
   * Per-thread number of cells already skipped in earlier hierarchy levels during getParticleImpl().
   */
  mutable std::vector<size_t> _prevCells;
  /**
   * Cache-line padding factor for _currentLevel and _prevCells to reduce false sharing between threads.
   */
  const size_t _cacheOffset;

  std::vector<std::unique_ptr<LinkedCells<Particle_T>>> _levels;
  std::vector<double> _maxCutoffPerLevel;
  /**
   * Indicates whether grids fit perfectly into one another, allowing for HGC08
   */
  bool _fittedGrids;
  /**
   * Gets the level of the hierarchical grid to which the particle should belong, e.g. for purposes of adding the
   * particle to the container. This is the level with the smallest cell size (excluding skin) that is greater than or
   * equal to that of the particle.
   * @param p the particle
   * @return the level of the hierarchical grid to which the particle should belong.
   */
  size_t getHierarchyLevel(const Particle_T &p) const {
    const double cutoff = p.getSize();
    for (size_t i = 0; i < _numLevels; ++i) {
      if (_maxCutoffPerLevel[i] >= cutoff) {
        return i;
      }
    }
    AutoPasLog(ERROR,
               "Size of Particle({}) is bigger than biggest allowed cutoff({}) of HierarchicalGrid, "
               "will result in wrong interaction calculation.",
               cutoff, _maxCutoffPerLevel.back());
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
      hGridTraversal->setLevels(&_levels, _maxCutoffPerLevel, _skin);
    } else {
      autopas::utils::ExceptionHandler::exception(
          "The selected traversal is not compatible with the HierarchicalGrid container. TraversalID: {}",
          traversal->getTraversalType());
    }
  }
};

}  // namespace autopas