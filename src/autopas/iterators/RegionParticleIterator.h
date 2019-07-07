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

namespace autopas::internal {
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
   * @param startRegion Lower corner of the region to iterate over.
   * @param endRegion Top corner of the region to iterate over.
   * @param indicesInRegion List of indices all threads will iterate over.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   */
  explicit RegionParticleIterator(std::vector<ParticleCell> *cont, std::array<double, 3> startRegion,
                                  std::array<double, 3> endRegion, std::vector<size_t> &indicesInRegion,
                                  internal::CellBorderAndFlagManager *flagManager = nullptr,
                                  IteratorBehavior behavior = haloAndOwned)
      : ParticleIterator<Particle, ParticleCell>(cont, flagManager, behavior),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(indicesInRegion),
        _currentRegionIndex(autopas_get_thread_num()) {
    this->_vectorOfCells = cont;
    if (_indicesInRegion.size() > (size_t)autopas_get_thread_num()) {
      this->_iteratorAcrossCells = cont->begin() + _indicesInRegion[autopas_get_thread_num()];
      this->_iteratorWithinOneCell = this->_iteratorAcrossCells->begin();
    } else {
      this->_iteratorAcrossCells = cont->end();
      return;
    }

    if (not this->isCellTypeBehaviorCorrect()) {
      this->next_non_empty_cell();
    }

    // The iterator might still be invalid (because the cell is empty or the owned-state of the particle is wrong), so
    // we check it here!
    if (ParticleIterator<Particle, ParticleCell>::isValid()) {
      // The iterator is valid, so there is a particle, now we need to check whether it's actually in the box!
      if (utils::notInBox(this->operator*().getR(), _startRegion, _endRegion)) {
        operator++();
      }
    } else if (this->_iteratorAcrossCells != cont->end()) {
      // the iterator is invalid + we still have particles, so we increment it!
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
      // make the iterator invalid!
      this->_iteratorAcrossCells = this->_vectorOfCells->end();
      return;
    }
    size_t iteratorInc = _indicesInRegion[_currentRegionIndex + stride] - _indicesInRegion[_currentRegionIndex];

    for (this->_iteratorAcrossCells += iteratorInc; this->getCurrentCellId() <= *(_indicesInRegion.end() - 1);
         this->_iteratorAcrossCells += iteratorInc) {
      _currentRegionIndex += stride;

      if (this->_iteratorAcrossCells < this->_vectorOfCells->end() and this->_iteratorAcrossCells->isNotEmpty() and
          this->isCellTypeBehaviorCorrect()) {
        this->_iteratorWithinOneCell = this->_iteratorAcrossCells->begin();
        break;
      }
      if (_currentRegionIndex + stride >= _indicesInRegion.size()) {
        // make the iterator invalid!
        this->_iteratorAcrossCells = this->_vectorOfCells->end();
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
}  // namespace autopas::internal