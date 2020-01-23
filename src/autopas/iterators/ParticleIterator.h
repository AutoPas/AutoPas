/**
 * @file ParticleIterator.h
 *
 * @authors tchievn, seckler
 * @date 17.01.2018
 */

#pragma once

#include <vector>

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/iterators/ParticleIteratorInterface.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/iterators/SingleCellIteratorWrapper.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas::internal {
/**
 * ParticleIterator class to access particles inside of a container.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator".
 * @tparam Particle Type of the particle that is accessed.
 * @tparam ParticleCell Type of the cell of the underlying data structure of the container.
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, class ParticleCell, bool modifiable>
class ParticleIterator : public ParticleIteratorInterfaceImpl<Particle, modifiable> {
  using CellVecType = std::conditional_t<modifiable, std::vector<ParticleCell>, const std::vector<ParticleCell>>;
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;
  using CellBorderAndFlagManagerType =
      std::conditional_t<modifiable, internal::CellBorderAndFlagManager, const internal::CellBorderAndFlagManager>;

 protected:
  /**
   * Simple constructor without major calculation overhead for internal use.
   * @note ATTENTION: This Iterator might be invalid after construction!
   *
   * @param cont Linear data vector of ParticleCells.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   */
  explicit ParticleIterator(CellVecType *cont, CellBorderAndFlagManagerType *flagManager, IteratorBehavior behavior)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell(cont->begin()->begin()),
        _flagManager(flagManager),
        _behavior(behavior) {}

 public:
  /**
   * Constructor of the ParticleIterator class.
   *
   * When started in a parallel block iterators are generated with offset
   * and iterate in a round robin fashion. If there are more threads than cells
   * surplus iterators are directly set to cont->end().
   *
   * @param cont Linear data vector of ParticleCells.
   * @param offset Number of cells to skip before starting to iterate.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   */
  explicit ParticleIterator(CellVecType *cont, size_t offset = 0, CellBorderAndFlagManagerType *flagManager = nullptr,
                            IteratorBehavior behavior = haloAndOwned)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell(cont->begin()->begin()),
        _flagManager(flagManager),
        _behavior(behavior) {
    auto myThreadId = autopas_get_thread_num();
    offset += myThreadId;
    if (offset < cont->size()) {
      _iteratorAcrossCells += offset;
      _iteratorWithinOneCell = _iteratorAcrossCells->begin();
    } else {
      _iteratorAcrossCells = cont->end();
      AutoPasLog(warn, "More threads than cells. No work left for thread {}!", myThreadId);
      return;
    }

    if (behavior != haloAndOwned and flagManager == nullptr) {
      AutoPasLog(error,
                 "Behavior is not haloAndOwned, but flagManager is "
                 "nullptr!");
      utils::ExceptionHandler::exception(
          "Behavior is not haloAndOwned, but flagManager is "
          "nullptr!");
    }

    if (not _iteratorWithinOneCell.isValid() or not isCellTypeBehaviorCorrect()) {
      ParticleIterator::next_non_empty_cell();
    }
    // The iterator might still be invalid (because the cell is empty or the owned-state of the particle is wrong), so
    // we check it here!
    if (not ParticleIterator::isValid() and this->_iteratorAcrossCells < cont->end()) {
      // the iterator is invalid + we still have particles/cells, so we increment it!
      ParticleIterator::operator++();
    }
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline ParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    // this iteration is done until a proper particle is found. i.e. one that satisfies the owning state.
    do {
      if (_iteratorWithinOneCell.isValid()) {
        ++_iteratorWithinOneCell;
      }

      // don't merge into if-else, _cell_iterator may become invalid after ++
      if (not _iteratorWithinOneCell.isValid()) {
        next_non_empty_cell();
      }
    } while (not isValid() && _iteratorAcrossCells < _vectorOfCells->end());
    return *this;
  }

  /**
   * @copydoc ParticleIteratorInterface::operator*()
   */
  inline ParticleType &operator*() const override { return _iteratorWithinOneCell.operator*(); }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  bool isValid() const override {
    return _vectorOfCells != nullptr and _iteratorAcrossCells < _vectorOfCells->end() and
           _iteratorWithinOneCell.isValid() and particleHasCorrectOwnedState();
  }

  ParticleIteratorInterfaceImpl<Particle, modifiable> *clone() const override {
    return new ParticleIterator<Particle, ParticleCell, modifiable>(*this);
  }

 protected:
  /**
   * @copydoc ParticleIteratorInterface::deleteCurrentParticle()
   */
  inline void deleteCurrentParticleImpl() override {
    if (_iteratorWithinOneCell.isValid()) {
      if constexpr (modifiable) {
        internal::deleteParticle(_iteratorWithinOneCell);
      } else {
        utils::ExceptionHandler::exception("Error: Trying to delete a particle through a const iterator.");
      }
    } else {
      utils::ExceptionHandler::exception("ParticleIterator: deleting invalid particle");
    }
  }

  /**
   * Moves the iterator to the next non-empty cell with respect to the stride.
   */
  virtual void next_non_empty_cell() {
    // find the next non-empty cell
    const int stride = autopas_get_num_threads();  // num threads
    for (_iteratorAcrossCells += stride; _iteratorAcrossCells < _vectorOfCells->end(); _iteratorAcrossCells += stride) {
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
  bool isCellTypeBehaviorCorrect() const {
    switch (_behavior) {
      case haloAndOwned:
        return true;
      case haloOnly:
        return _flagManager->cellCanContainHaloParticles(_iteratorAcrossCells - _vectorOfCells->begin());
      case ownedOnly:
        return _flagManager->cellCanContainOwnedParticles(_iteratorAcrossCells - _vectorOfCells->begin());
      default:
        utils::ExceptionHandler::exception("unknown iterator behavior");
        return false;
    }
  }

  /**
   * Indicates whether the particle has the correct owned state.
   * @return
   */
  bool particleHasCorrectOwnedState() const {
    switch (_behavior) {
      case haloAndOwned:
        return true;
      case haloOnly:
        return not _iteratorWithinOneCell->isOwned();
      case ownedOnly:
        return _iteratorWithinOneCell->isOwned();
      default:
        utils::ExceptionHandler::exception("unknown iterator behavior");
        return false;
    }
  }

  /**
   * Get the 1D index of the cell the iterator currently is in.
   * @return Index of current cell.
   */
  size_t getCurrentCellId() { return _iteratorAcrossCells - _vectorOfCells->begin(); }

  /**
   * Pointer to the cell vector.
   */
  CellVecType *_vectorOfCells;

  /**
   * Iterator for traversing the cell vector.
   */
  std::conditional_t<modifiable, typename CellVecType::iterator, typename CellVecType::const_iterator>
      _iteratorAcrossCells;

  /**
   * Particle iterator for a single cell.
   */
  SingleCellIteratorWrapper<Particle, modifiable> _iteratorWithinOneCell;

  /**
   * Manager providing info if cell is in halo.
   */
  CellBorderAndFlagManagerType *_flagManager;

  /**
   * The behavior of the iterator.
   */
  const IteratorBehavior _behavior;
};
}  // namespace autopas::internal
