/**
 * @file LinkedCells.h
 * @date 05.04.2020
 * @author lunaticcoding
 */

#pragma once

#include <autopas/utils/Timer.h>

#include "autopas/cells/ReferenceParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
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
#include <algorithm>

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

  enum SortingAlgorithm {
    none, stdsort, quicksort, radixsort, countingsort
  };
  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle;
  /**
   *  Type of the ParticleCell.
   */
  using ReferenceCell = ReferenceParticleCell<Particle>;

  SortingAlgorithm sorting;
  int sort_loop_iterations;
  std::array<double, 3> boxmin_, boxmax_;
  autopas::utils::Timer sorting_timer;
  /**
   * Constructor of the LinkedCells class
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   * By default all applicable traversals are allowed.
   */
  LinkedCellsReferences(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                        const double skin, const double cellSizeFactor = 1.0,
                        LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell, SortingAlgorithm _sorting = quicksort, int _sort_loop_iterations = 100)
      : CellBasedParticleContainer<ReferenceCell>(boxMin, boxMax, cutoff, skin),
        _cellBlock(this->_cells, boxMin, boxMax, cutoff + skin, cellSizeFactor),
        _loadEstimator(loadEstimator) {
    sorting = _sorting;
    boxmin_ = boxMin;
    boxmax_ = boxMax;
    sort_loop_iterations = _sort_loop_iterations;
  }

  double keyof(Particle x){
//    return x.getR()[0];
    return static_cast<double>(_cellBlock.get1DIndexOfPosition(x.getR()));
  }
  bool compare(Particle a, Particle b){
    return keyof(a) < keyof(b);
  }
  enum QuicksortPivotStrategy {
    first, middle, med3, last
  };

  void sort(){
    sorting_timer.start();
    switch(sorting){
      case stdsort: {
        sort_std();
        break;
      }
      case none:
        break;
      case countingsort: {
        throw "Counting sort not implemented yet.";
        break;
      }
      case quicksort: {
        sort_quick();
        break;
      }
      case radixsort: {
        throw "Radix sort not implemented yet.";
        break;
      }
    }
    // updateContainer();
    sorting_timer.stop();
    std::cout << "Sorting took "<<sorting_timer.getTotalTime()<<" in total so far."<<std::endl;
  }
  void sort_std(){
    LinkedCellsReferences* lcr = this;
    std::sort(_particleList.begin(), _particleList.end(),
              [&lcr](const Particle a, const Particle b) -> bool { return lcr->compare(a,b) ; });
  }
  int64_t getMantissa(double d){
    return reinterpret_cast<int64_t&>(d) & (((int64_t)1 << 53) - 1);
  }
  int getExp(double d){
    int exp;
    frexp(d, &exp);
    return exp;
  }
  void sort_radix_simple(){
    // currently only an extremely reduced version to only approximate the effect of sorting.
    // (also currently not implemented because of in-place demand)
    throw "Radix sort not yet implemented.";
  }
  void sort_radix(){
    // for cache efficiency, it should be enough to only sort some bits.
    // bit_accuracy is the number of MSDs that are used for sorting. This can be changed empirically!:
    int bit_accuracy = 10;
    // ...also possible to set it depending on the approx. bit differences of positions inside different cells
    int maxexp = getExp(boxmax_[0]);
    int minexp = getExp(boxmin_[0]);
    int n = 0;
    for(int filter = 0x8000; n < 16 && (maxexp & filter) == (minexp & filter); filter >>= 1, n++);
    // first key: exponent, first position n where exponents differ in [0,16] (16 if exponent is equal)
    // (later since exponent contains MSDs)
    // (16-n) are the bits we need to handle with the first key
    // the remaining (bit_accuracy - (16-n)) bits will be handled with the second key
    // second key: mantissa, only need to inspect (bit_accuracy - (16-n)) bits (if <= 0, nothing has to be done)
    int second_n = bit_accuracy - (16 - n);
    if(second_n > 0) {
      // mantissa is 53 bit
      int64_t mantissa_filter = ((int64_t)1) << 52;
      // special case: min and max exponents are equal. Then:
      if(second_n == bit_accuracy) { // search first bit difference in mantissa
        int64_t maxmant = getMantissa(boxmax_[0]);
        int64_t minmant = getMantissa(boxmin_[0]);
        int n2 = 0;
        for (; n2 < 53 && (maxmant & mantissa_filter) == (minmant & mantissa_filter); mantissa_filter >>= 1, n2++);
        int still_remaining = (53 - n2); // n2=53 iff mantissas are equal (iff boxmin[0]=boxmax[0])
        if(still_remaining < bit_accuracy){
          // accuracy is less than bit_accuracy -> fewer digits to sort!
          second_n -= (bit_accuracy - still_remaining);
        }
      }
      // radix sort stage 1 - mantissa
      sort_radix_mantissa(second_n, mantissa_filter);

    }

  }
  void swap(size_t i, size_t j){
//       ParticleType temp = _particleList[i];
//       _particleList[i] = _particleList[j];
//       _particleList[j] = temp;
      std::swap((_particleList[i]), (_particleList[j]));
  }
  bool get_nth_bit_mantissa(double d, int n){
    // true iff 1
    int64_t mantissa = getMantissa(d);
    int64_t mask = ((int64_t)1) << (52 - n);
    return (mantissa & mask) == mask;
  }
  bool get_nth_bit_exponent_16bit(double d, int n){
    int exp = getExp(d);
    int mask = 1 << (15 - n);
    return (exp & mask) == mask;
  }
  void sort_radix_mantissa(int digits, int64_t mantissa_filter = ((int64_t)1) << 52){
    throw "Radix Sort: binary digit sorting -> mantissa sorting stage not yet implemented.";
    for(; digits > 0; digits--, mantissa_filter >>= 1){
      int i = 0;
      int j = _particleList.totalSize() - 1;
      while (i < j){
        //while (i < j && _particleList[i]);
      }
    }
  }
  void sort_radix_exp(){
    throw "Radix Sort: binary digit sorting -> exponent sorting stage not yet implemented.";
  }
  void sort_quick(QuicksortPivotStrategy pivot = last, size_t left = 0, size_t right = 0, bool top_level = true){
    if(_particleList.totalSize() == 0) return;
    if(top_level)
      right = _particleList.totalSize() - 1;
    if (right <= left) return;
    switch(pivot){
      case last: {
        break;
      }
      case first: {
        swap(left, right);
        break;
      }
      case middle: {
        size_t m = (left+right)/2;
        swap(m, right);
        break;
      }
      case med3: {
        throw "3-Median not implemented yet.";
        break;
      }
    }
    double pivotkey = keyof(_particleList[right]);
    // Hoare partitioning
    size_t i = left-1, j = right;
    do {
      do { i++; } while (keyof(_particleList[i]) < pivotkey);
      do {
        if(j == 0) break;
        j--;
      } while (j >= left && keyof(_particleList[j]) > pivotkey);
      if(i < j) swap(i,j);
    } while (i < j);
    this->swap(right, i);
    if(i != 0) //i == 0 -> size_t underflow. Solution: skip since nothing has to be done for index interval [0,0]
      sort_quick(pivot, left, i-1, false);
    sort_quick(pivot, i+1, right, false);
  }
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
    auto cells = _cellBlock.getNearbyHaloCells(pCopy.getR(), this->getSkin());
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
//    if(sort_counter % sort_loop_iterations == 0){
//      sort();
//    }
//    sort_counter++;
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

  std::vector<ParticleType> updateContainer() override {
    this->deleteHaloParticles();

    // if sorting is none, the following will have no impact
    sort();

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
        if (not this->getCells()[cellId].isNotEmpty()) continue;

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

  ParticleIteratorWrapper<ParticleType, true> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                const std::array<double, 3> &higherCorner,
                                                                IteratorBehavior behavior) override {
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
        new internal::RegionParticleIterator<ParticleType, ReferenceCell, true>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior, nullptr));
  }

  ParticleIteratorWrapper<ParticleType, false> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                                 const std::array<double, 3> &higherCorner,
                                                                 IteratorBehavior behavior) const override {
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
        new internal::RegionParticleIterator<ParticleType, ReferenceCell, false>(
            &this->_cells, lowerCorner, higherCorner, cellsOfInterest, &_cellBlock, behavior, nullptr));
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

  int sort_counter = 0;
};

}  // namespace autopas
