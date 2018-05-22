/*
 * ParticleIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_

#include <utils/ExceptionHandler.h>
#include <vector>
#include "SingleCellIterator.h"

namespace autopas {

/**
 * Enum to specify the behavior of an iterator.
 */
enum IteratorBehavior {
  haloOnly,     /// iterate only over halo
  ownedOnly,    /// iterate only over inner cells
  haloAndOwned  /// iterate over both halo and inner cells
};

/**
 * ParticleIterator class to access particles inside of a container.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator"
 * @tparam Particle type of the particle that is accessed
 * @tparam ParticleCell type of the cell of the underlying data structure of the
 * container
 */
template <class Particle, class ParticleCell>
class ParticleIterator {
 public:
  /**
   * Constructor of the ParticleIterator class.
   *
   * @param cont linear data vector of ParticleCells
   * @param flagManager the CellBorderAndFlagManager that shall be used to
   * query the cell types. Can be nullptr if the behavior is haloAndOwned
   * @param behavior the IteratorBehavior that specifies which type of cells
   * shall be iterated through.
   */
  explicit ParticleIterator(std::vector<ParticleCell>* cont,
                            CellBorderAndFlagManager* flagManager = nullptr,
                            IteratorBehavior behavior = haloAndOwned)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell(),
        _flagManager(flagManager),
        _behavior(behavior) {
    if (behavior != haloAndOwned and flagManager == nullptr) {
      AutoPasLogger->error(
          "ParticleIterator: behavior is not haloAndOwned, but flagManager is "
          "nullptr!");
      utils::ExceptionHandler::exception(
          "ParticleIterator: behavior is not haloAndOwned, but flagManager is "
          "nullptr!");
    }
    if (_iteratorAcrossCells < cont->end()) {
      _iteratorWithinOneCell = _iteratorAcrossCells->begin();
      if (not isCellTypeBehaviorCorrect() or
          not _iteratorWithinOneCell.isValid()) {
        next_non_empty_cell();
      }
    }
  }

  /**
   * Increment operator.
   * Used to jump to the next particle
   * @return next particle, usually ignored
   */
  inline ParticleIterator<Particle, ParticleCell>& operator++() {
    if (_iteratorWithinOneCell.isValid()) {
      ++_iteratorWithinOneCell;
    }

    // don't merge into if-else, _cell_iterator may becoeme invalid after ++

    if (not _iteratorWithinOneCell.isValid()) {
      next_non_empty_cell();
    }
    return *this;
  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  Particle& operator*() { return _iteratorWithinOneCell.operator*(); }

  /**
   * access particle using iterator->
   *
   * this is the member of pointer operator
   * @return current particle
   */
  Particle* operator->() { return &(this->operator*()); }  //

  /**
   * Deletes the current particle
   */
  void deleteCurrentParticle() {
    if (_iteratorWithinOneCell.isValid()) {
      _iteratorWithinOneCell.deleteCurrentParticle();
    } else {
      utils::ExceptionHandler::exception(
          "ParticleIterator: deleting invalid particle");
    }
  }

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  bool isValid() {
    return _vectorOfCells != nullptr and
           _iteratorAcrossCells < _vectorOfCells->end() and
           _iteratorWithinOneCell.isValid();
  }

 protected:
  /**
   * the next non-empty cell is selected
   */
  void next_non_empty_cell() {
    // find the next non-empty cell
    const int stride = 1;  // num threads
    for (_iteratorAcrossCells += stride;
         _iteratorAcrossCells < _vectorOfCells->end();
         _iteratorAcrossCells += stride) {
      if (_iteratorAcrossCells->isNotEmpty() and isCellTypeBehaviorCorrect()) {
        _iteratorWithinOneCell = _iteratorAcrossCells->begin();
        break;
      }
    }
  }

  /**
   * checks if a cell has the correct cell type according to the behavior
   * @return true iff the cell type is proper according to the behavior
   */
  bool isCellTypeBehaviorCorrect() {
    switch (_behavior) {
      case haloAndOwned:
        return true;
      case haloOnly:
        return _flagManager->isHaloCell(_iteratorAcrossCells -
                                        _vectorOfCells->begin());
      case ownedOnly:
        return _flagManager->isOwningCell(_iteratorAcrossCells -
                                          _vectorOfCells->begin());
      default:
        return false;
    }
  }

 private:
  std::vector<ParticleCell>* _vectorOfCells;
  typename std::vector<ParticleCell>::iterator _iteratorAcrossCells;
  typename ParticleCell::iterator _iteratorWithinOneCell;
  CellBorderAndFlagManager* _flagManager;
  IteratorBehavior _behavior;
  // SingleCellIterator<Particle, ParticleCell> _iteratorWithinOneCell;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_ */
