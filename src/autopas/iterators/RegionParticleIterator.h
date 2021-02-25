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
 * @tparam Particle Particle type over which the iterator iterates
 * @tparam ParticleCell Cell type over which the iterator iterates
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, class ParticleCell, bool modifiable>
class RegionParticleIterator : public ParticleIterator<Particle, ParticleCell, modifiable> {
  using CellVecType = std::conditional_t<modifiable, std::vector<ParticleCell>, const std::vector<ParticleCell>>;
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;
  using CellBorderAndFlagManagerType =
      std::conditional_t<modifiable, internal::CellBorderAndFlagManager, const internal::CellBorderAndFlagManager>;
  using ParticleVecType =
      std::conditional_t<modifiable, std::vector<std::vector<Particle>>, const std::vector<std::vector<Particle>>>;

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
   * @param additionalParticleVectorToIterate Additional Particle Vector to iterate over.
   */
  explicit RegionParticleIterator(CellVecType *cont, std::array<double, 3> startRegion, std::array<double, 3> endRegion,
                                  std::vector<size_t> &indicesInRegion,
                                  const CellBorderAndFlagManagerType *flagManager = nullptr,
                                  IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                  ParticleVecType *additionalParticleVectorToIterate = nullptr)
      : ParticleIterator<Particle, ParticleCell, modifiable>(cont, flagManager, behavior,
                                                             additionalParticleVectorToIterate),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(indicesInRegion) {
    auto myThreadId = autopas_get_thread_num();
    _currentRegionIndex = myThreadId;
    if (additionalParticleVectorToIterate and myThreadId == autopas_get_num_threads() - 1) {
      // we want to iterate with the last thread over the additional particle vector.
      this->_additionalParticleVectorToIterateState =
          decltype(this->_additionalParticleVectorToIterateState)::notStarted;
    }

    this->_vectorOfCells = cont;
    if (_indicesInRegion.size() > static_cast<size_t>(myThreadId)) {
      this->_iteratorAcrossCells = cont->begin() + _indicesInRegion[myThreadId];
      this->_iteratorWithinOneCell = this->_iteratorAcrossCells->begin();
    } else if (this->_additionalParticleVectorToIterateState ==
               decltype(this->_additionalParticleVectorToIterateState)::notStarted) {
      this->_additionalParticleVectorToIterateState =
          decltype(this->_additionalParticleVectorToIterateState)::iterating;
    } else {
      this->_iteratorAcrossCells = cont->end();
      return;
    }

    if (behavior != IteratorBehavior::haloAndOwned and flagManager == nullptr) {
      AutoPasLog(error,
                 "Behavior is not haloAndOwned, but flagManager is "
                 "nullptr!");
      utils::ExceptionHandler::exception(
          "Behavior is not haloAndOwned, but flagManager is "
          "nullptr!");
    }

    if (this->_additionalParticleVectorToIterateState !=
        decltype(this->_additionalParticleVectorToIterateState)::iterating) {
      if (not this->isCellTypeBehaviorCorrect()) {
        this->next_non_empty_cell();
      }

      // The iterator might still be invalid (because the cell is empty or the owned-state of the particle is wrong), so
      // we check it here!
      if (ParticleIterator<Particle, ParticleCell, modifiable>::isValid()) {
        // The iterator is valid, so there is a particle, now we need to check whether it's actually in the box!
        if (utils::notInBox(this->operator*().getR(), _startRegion, _endRegion)) {
          operator++();
        }
      } else if (this->_iteratorAcrossCells != cont->end()) {
        // the iterator is invalid + we still have particles, so we increment it!
        operator++();
      }
    } else {
      if (not isValid()) {
        operator++();
      }
    }
  }

  inline RegionParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    do {
      ParticleIterator<Particle, ParticleCell, modifiable>::operator++();
    } while (ParticleIterator<Particle, ParticleCell, modifiable>::isValid() and
             utils::notInBox(this->operator*().getR(), _startRegion, _endRegion) and
             this->getCurrentCellId() <= *(_indicesInRegion.end() - 1));
    return *this;
  }

  bool isValid() const override {
    return ParticleIterator<Particle, ParticleCell, modifiable>::isValid() &&
           utils::inBox(this->operator*().getR(), _startRegion, _endRegion);
  }

  inline ParticleIteratorInterfaceImpl<Particle, modifiable> *clone() const override {
    return new RegionParticleIterator<Particle, ParticleCell, modifiable>(*this);
  }

 private:
  /**
   * @copydoc ParticleIterator::next_non_empty_cell()
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

  const std::array<double, 3> _startRegion;
  const std::array<double, 3> _endRegion;
  const std::vector<size_t> _indicesInRegion;
  size_t _currentRegionIndex;
};
}  // namespace autopas::internal
