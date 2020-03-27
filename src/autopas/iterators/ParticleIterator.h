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
  using ParticleVecType = std::conditional_t<modifiable, std::vector<Particle>, const std::vector<Particle>>;

 protected:
  /**
   * Simple constructor without major calculation overhead for internal use.
   * @note ATTENTION: This Iterator might be invalid after construction!
   *
   * @param cont Linear data vector of ParticleCells.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   * @param additionalParticleVectorToIterate Additional Particle Vector to iterate over.
   */
  ParticleIterator(CellVecType *cont, CellBorderAndFlagManagerType *flagManager, IteratorBehavior behavior,
                   ParticleVecType *additionalParticleVectorToIterate)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell(cont->begin()->begin()),
        _flagManager(flagManager),
        _behavior(behavior),
        _additionalParticleVector(additionalParticleVectorToIterate) {}

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
   * @param additionalParticleVectorToIterate Additional Particle Vector to iterate over.
   */
  explicit ParticleIterator(CellVecType *cont, size_t offset = 0, CellBorderAndFlagManagerType *flagManager = nullptr,
                            IteratorBehavior behavior = haloAndOwned,
                            ParticleVecType *additionalParticleVectorToIterate = nullptr)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell(cont->begin()->begin()),
        _flagManager(flagManager),
        _behavior(behavior),
        _additionalParticleVector{additionalParticleVectorToIterate} {
    auto myThreadId = autopas_get_thread_num();
    offset += myThreadId;
    if (additionalParticleVectorToIterate and myThreadId == autopas_get_num_threads() - 1) {
      // we want to iterate with the last thread over the additional particle vector.
      _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::notStarted;
    }

    if (offset < cont->size()) {
      _iteratorAcrossCells += offset;
      _iteratorWithinOneCell = _iteratorAcrossCells->begin();
    } else if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::notStarted) {
      _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::iterating;
    } else {
      _iteratorAcrossCells = cont->end();
      AutoPasLog(trace, "More threads than cells. No work left for thread {}!", myThreadId);
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

    if (_additionalParticleVectorToIterateState != AdditionalParticleVectorToIterateState::iterating) {
      if (not _iteratorWithinOneCell.isValid() or not isCellTypeBehaviorCorrect()) {
        // This ensures that the current cell is NOT empty or (that the end of the cells is selected and
        // ParticleIterator::isValid is false).
        /// @todo: this might not actually be needed -> check!
        ParticleIterator::next_non_empty_cell();
      }
      // The iterator might still be invalid (because the cell is empty or the owned-state of the particle is wrong), so
      // we check it here!
      if (not ParticleIterator::isValid() and
          (this->_iteratorAcrossCells < cont->end() or
           _additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::notStarted)) {
        // the iterator is invalid + we still have particles/cells, so we increment it!
        ParticleIterator::operator++();
      }
    } else {
      if (not ParticleIterator::isValid()) {
        // The iterator is invalid + we still have particles/cells, so we increment it!
        ParticleIterator::operator++();
      }
    }
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline ParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    if (_additionalParticleVectorToIterateState != AdditionalParticleVectorToIterateState::iterating) {
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
      if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::notStarted &&
          _iteratorAcrossCells >= _vectorOfCells->end()) {
        _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::iterating;
        if (isValid()) {
          // if the first particle in the additional Vector is valid, we return here.
          return *this;
        }  // otherwise we run into the loop below and ++ further.
      }
    }
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      // Simply increase the counter, as long as we aren't valid.
      do {
        ++_additionalParticleVectorPosition;
      } while (not isValid() && _additionalParticleVectorPosition < _additionalParticleVector->size());
    }

    return *this;
  }

  /**
   * @copydoc ParticleIteratorInterface::operator*()
   */
  inline ParticleType &operator*() const override {
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      return (*_additionalParticleVector)[_additionalParticleVectorPosition];
    }
    return _iteratorWithinOneCell.operator*();
  }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  bool isValid() const override {
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      return _additionalParticleVectorPosition < _additionalParticleVector->size() and particleHasCorrectOwnedState();
    }
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
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      if constexpr (modifiable) {
        // for the additionalParticleVector, we simply swap the last particle with the current particle and remove the
        // particle that is now at the back.
        std::swap((*_additionalParticleVector)[_additionalParticleVectorPosition], _additionalParticleVector->back());
        _additionalParticleVector->pop_back();
        --_additionalParticleVectorPosition;
      } else {
        utils::ExceptionHandler::exception("Error: Trying to delete a particle through a const iterator.");
      }
    } else if (_iteratorWithinOneCell.isValid()) {
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
        if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
          return not(*_additionalParticleVector)[_additionalParticleVectorPosition].isOwned();
        } else {
          return not _iteratorWithinOneCell->isOwned();
        }
      case ownedOnly:
        if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
          return (*_additionalParticleVector)[_additionalParticleVectorPosition].isOwned();
        } else {
          return _iteratorWithinOneCell->isOwned();
        }
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

  /**
   * Enum class that specifies the iterating state of the additional particle vector.
   */
  enum class AdditionalParticleVectorToIterateState { ignore, notStarted, iterating };

  /**
   * Specifies if this iterator should also iterate over the additional array.
   */
  AdditionalParticleVectorToIterateState _additionalParticleVectorToIterateState{
      AdditionalParticleVectorToIterateState::ignore};

  /**
   * Pointer to an additional Particle Vector this ParticleIterator will iterate over.
   */
  ParticleVecType *_additionalParticleVector{nullptr};

  /**
   * Index for ParticleVecType.
   */
  size_t _additionalParticleVectorPosition{0ul};
};
}  // namespace autopas::internal
