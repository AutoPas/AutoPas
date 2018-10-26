/**
 * @file RegionParticleIterator.h specifies the RegionParticleIterator class
 * @author seckler
 * @date 03.04.2018
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/inBox.h"

namespace autopas {
namespace internal {
/**
 * RegionParticleIterator to iterate over all particles within a specific region
 * @todo optimize the region particle iterater. Currently we iterate over all
 * particles
 * @tparam Particle Particle type over which the iterator iterates
 * @tparam ParticleCell Cell type over which the iterator iterates
 */
template <class Particle, class ParticleCell>
class RegionParticleIterator : public ParticleIterator<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the RegionParticleIterator.
   *
   * @param cont Container of particle cells.
   * @param cellsPerDimension Total number of cells along each dimension.
   * @param startRegion Lower corner of the region to iterate over.
   * @param endRegion Top corner of the region to iterate over.
   * @param startIndex Index of the first cell in the region of interest.
   * @param endIndex Index of the last cell in the region of interest.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   */
  explicit RegionParticleIterator(std::vector<ParticleCell> *cont, std::array<double, 3> startRegion,
                                  std::array<double, 3> endRegion, std::vector<size_t> &indicesInRegion,
                                  CellBorderAndFlagManager *flagManager = nullptr,
                                  IteratorBehavior behavior = haloAndOwned)
      : ParticleIterator<Particle, ParticleCell>(cont, flagManager, behavior),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(indicesInRegion),
        _currentRegionIndex(autopas_get_thread_num()) {
    this->_vectorOfCells = cont;
    if (_indicesInRegion.size() >= (size_t)autopas_get_num_threads()) {
      this->_iteratorAcrossCells = cont->begin() + _indicesInRegion[autopas_get_thread_num()];
      this->_iteratorWithinOneCell = this->_iteratorAcrossCells->begin();
    } else {
      this->_iteratorAcrossCells = cont->end();
    }

    this->_flagManager = flagManager;
    this->_behavior = behavior;
    // ParticleIterator's constructor will initialize the Iterator, such that it
    // points to the first particle if one is found, otherwise the pointer is
    // not valid
    if (ParticleIterator<Particle, ParticleCell>::isValid()) {  // if there is NO particle, we can not dereference
                                                                // it, so we need a check.
      if (utils::notInBox(this->operator*().getR(), _startRegion, _endRegion)) {
        operator++();
      }
    } else if (this->_iteratorAcrossCells != cont->end()) {
      operator++();
    }
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++
   */
  inline RegionParticleIterator<Particle, ParticleCell> &operator++() override {
    do {
      ParticleIterator<Particle, ParticleCell>::operator++();
    } while (ParticleIterator<Particle, ParticleCell>::isValid() &&
             utils::notInBox(this->operator*().getR(), _startRegion, _endRegion) &&
             this->getCurrentCellId() <= *(_indicesInRegion.end() - 1));
    return *this;
  }

  bool isValid() const override {
    return ParticleIterator<Particle, ParticleCell>::isValid() &&
           utils::inBox(this->operator*().getR(), _startRegion, _endRegion);
  }

  // @todo add test of clone
  inline ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new RegionParticleIterator<Particle, ParticleCell>(*this);
  }

 private:
  /**
   * @copydoc
   * @note overrides call in ParticleIterator<Particle, ParticleCell>::operator++()
   */
  void next_non_empty_cell() override {
    // find the next non-empty cell
    const int stride = autopas_get_num_threads();  // num threads
    if (_currentRegionIndex + stride >= _indicesInRegion.size()) {
      return;
    }
    size_t iteratorInc = _indicesInRegion[_currentRegionIndex + stride] - _indicesInRegion[_currentRegionIndex];

    for (this->_iteratorAcrossCells += iteratorInc; this->getCurrentCellId() <= *(_indicesInRegion.end() - 1);
         this->_iteratorAcrossCells += iteratorInc) {
      //#pragma omp critical
      //      std::cout << "Thread [" << autopas_get_thread_num() << "] currentCellIndex= " <<
      //      (this->getCurrentCellId())
      //      << " currentRegionIndex= " << _currentRegionIndex
      //      << " iteratorInc= " << iteratorInc << std::endl;
      _currentRegionIndex += stride;

      if (this->_iteratorAcrossCells < this->_vectorOfCells->end() and this->_iteratorAcrossCells->isNotEmpty() and
          this->isCellTypeBehaviorCorrect()) {
        this->_iteratorWithinOneCell = this->_iteratorAcrossCells->begin();
        break;
      }
      if (_currentRegionIndex + stride >= _indicesInRegion.size()) {
        break;
      }
      iteratorInc = _indicesInRegion[_currentRegionIndex + stride] - _indicesInRegion[_currentRegionIndex];
    }
  }

  std::array<double, 3> _startRegion;
  std::array<double, 3> _endRegion;
  std::vector<size_t> _indicesInRegion;
  size_t _currentRegionIndex;
};
}  // namespace internal
}  // namespace autopas