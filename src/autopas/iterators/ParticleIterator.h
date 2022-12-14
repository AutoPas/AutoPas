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
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas::internal {
/**
 * ParticleIterator class to access particles inside of a cell based container.
 * The particles can be accessed using "iterator->" or "*iterator". The iterator is incremented using the ++operator.
 *
 * Incrementing this iterator works by jumping through the underlying particle cells from particle to particle
 * incrementing until it finds a particle that matches the iteration criteria (e.g. it continues to loop until it finds
 * a halo particle if _behavior==halo). Additionally, there might be additional particle vectors that are not part
 * of the cell structure and contain arbitrary unsorted particles. These buffers result from containers that
 * have buffers they only sometimes merge with the actual data structure. Typically there is one buffer per OpenMP
 * thread, however, there is no guarantee that the iterator region has the same number of threads so the number of
 * buffers is considered arbitrary.
 *
 * When instantiated in a parallel region one iterator per thread is spawned. All iterators start at the cell of their
 * thread ID and then jump cells by the number of threads. When a thread is done with its particles it continues with
 * the additional buffers. Again starting at the buffer with its thread ID and jumping by the number of threads.
 *
 * @tparam Particle Type of the particle that is accessed.
 * @tparam ParticleCell Type of the cell of the underlying data structure of the container.
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, class ParticleCell, bool modifiable>
class ParticleIterator : public ParticleIteratorInterfaceImpl<Particle, modifiable> {
  using CellVecType = std::conditional_t<modifiable, std::vector<ParticleCell>, const std::vector<ParticleCell>>;
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;
  using CellBorderAndFlagManagerType = const internal::CellBorderAndFlagManager;
  using ParticleVecTypeInput =
      std::conditional_t<modifiable, std::vector<std::vector<Particle>>, const std::vector<std::vector<Particle>>>;
  using ParticleVecTypeInternal =
      std::conditional_t<modifiable, std::vector<std::vector<Particle> *>, std::vector<std::vector<Particle> const *>>;

 protected:
  /**
   * Simple constructor without major calculation overhead for internal use.
   * @note ATTENTION: This Iterator might be invalid after construction!
   *
   * @param cont Linear data vector of ParticleCells.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is ownedOrHalo.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   * @param additionalVectorsToIterate Thread buffers of additional Particle vector to iterate over.
   */
  ParticleIterator(CellVecType *cont, const CellBorderAndFlagManagerType *flagManager, IteratorBehavior behavior,
                   ParticleVecTypeInput *additionalVectorsToIterate)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell(cont->begin()->begin()),
        _flagManager(flagManager),
        _behavior(behavior),
        _additionalParticleVectorToIterateState(additionalVectorsToIterate
                                                    ? AdditionalParticleVectorToIterateState::notStarted
                                                    : AdditionalParticleVectorToIterateState::ignore),
        _additionalVectorIndex((behavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num()),
        _additionalVectorPosition(0) {
    if (additionalVectorsToIterate) {
      // heuristic: Owned and Halo vectors for each currently active thread
      _additionalVectors.reserve(autopas_get_num_threads() * 2);
      for (auto &&vec : *additionalVectorsToIterate) {
        _additionalVectors.push_back(&vec);
      }
    }
  }

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
   * Can be nullptr if the behavior is ownedOrHalo.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   * @param additionalVectorsToIterate Additional Particle Vector to iterate over.
   */
  explicit ParticleIterator(CellVecType *cont, size_t offset = 0, CellBorderAndFlagManagerType *flagManager = nullptr,
                            IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
                            ParticleVecTypeInput *additionalVectorsToIterate = nullptr)
      : ParticleIterator(cont, flagManager, behavior, additionalVectorsToIterate) {
    if (not(_behavior & IteratorBehavior::forceSequential)) {
      offset += autopas_get_thread_num();
    }

    if (offset < cont->size()) {
      _iteratorAcrossCells += offset;
      _iteratorWithinOneCell = _iteratorAcrossCells->begin();
    } else if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::notStarted) {
      _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::iterating;
    } else {
      _iteratorAcrossCells = cont->end();
      AutoPasLog(trace, "More threads than cells. No work left for thread {}!", autopas_get_thread_num());
      return;
    }

    if (not(behavior & IteratorBehavior::ownedOrHalo) and flagManager == nullptr) {
      utils::ExceptionHandler::exception(
          "ParticleIterator::ParticleIterator() Behavior is not ownedOrHalo, but flagManager is nullptr!");
    }

    if (_additionalParticleVectorToIterateState != AdditionalParticleVectorToIterateState::iterating) {
      if (not _iteratorWithinOneCell.isValid() or not isCellTypeBehaviorCorrect()) {
        // This ensures that the first cell has a correct cell type behavior. This couldn't be done using operator++.
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
      } while (not isValid() and _iteratorAcrossCells < _vectorOfCells->end());
      if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::notStarted and
          _iteratorAcrossCells >= _vectorOfCells->end() and _additionalVectors.size() > _additionalVectorIndex) {
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
        if (_additionalVectorIndex >= _additionalVectors.size()) {
          // return if we get to be invalid!
          _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::ignore;
          return *this;
        }
        ++_additionalVectorPosition;
        // if we reach the end of this buffer jump to the next
        if (_additionalVectorPosition >= _additionalVectors[_additionalVectorIndex]->size()) {
          _additionalVectorIndex += (_behavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();
          _additionalVectorPosition = 0;
        }
        // continue looking for valid iterator positions, while we aren't valid.
      } while (not isValid());
    }
    return *this;
  }

  /**
   * @copydoc ParticleIteratorInterface::operator*()
   */
  inline ParticleType &operator*() const override {
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      return (*_additionalVectors[_additionalVectorIndex])[_additionalVectorPosition];
    }
    return _iteratorWithinOneCell.operator*();
  }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  [[nodiscard]] bool isValid() const override {
    // case: currently iterating extra vectors
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      // true if: pointing at a valid particle within any of the vectors
      return _additionalVectorIndex < _additionalVectors.size() and
             _additionalVectorPosition < _additionalVectors[_additionalVectorIndex]->size() and
             particleHasCorrectOwnershipState();
    }
    // true if: iterator across and within cells are good
    return _vectorOfCells != nullptr and _iteratorAcrossCells < _vectorOfCells->end() and
           _iteratorWithinOneCell.isValid() and particleHasCorrectOwnershipState();
  }

  ParticleIteratorInterfaceImpl<Particle, modifiable> *clone() const override {
    return new ParticleIterator<Particle, ParticleCell, modifiable>(*this);
  }

  void addAdditionalVectors(
      std::conditional_t<modifiable, std::vector<std::vector<Particle>> &, const std::vector<std::vector<Particle>> &>
          additionalVectors) override {
    _additionalVectors.reserve(_additionalVectors.size() + additionalVectors.size());
    for (auto &vec : additionalVectors) {
      if (not vec.empty()) {
        _additionalVectors.push_back(&vec);
      }
    }
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::ignore) {
      _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::notStarted;
    }
    if (not isValid()) {
      // In case the iterator wasn't valid before, it might be now with the new data. Hence, increment it to find out.
      operator++();
    }
  }

  void addAdditionalVector(std::conditional_t<modifiable, std::vector<Particle> &, const std::vector<Particle> &>
                               additionalVector) override {
    _additionalVectors.push_back(&additionalVector);
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::ignore) {
      _additionalParticleVectorToIterateState = AdditionalParticleVectorToIterateState::notStarted;
    }
    if (not isValid()) {
      // In case the iterator wasn't valid before, it might be now with the new data. Hence, increment it to find out.
      operator++();
    }
  }

 protected:
  /**
   * @copydoc ParticleIteratorInterface::deleteCurrentParticleImpl()
   */
  inline void deleteCurrentParticleImpl() override {
    if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      if constexpr (modifiable) {
        auto &currentAdditionalVector = *_additionalVectors[_additionalVectorIndex];
        // for the additionalParticleVector, we simply swap the last particle with the current particle and remove the
        // particle that is now at the back.
        std::swap(currentAdditionalVector[_additionalVectorPosition], currentAdditionalVector.back());
        currentAdditionalVector.pop_back();
        --_additionalVectorPosition;
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
    const int stride = (_behavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();
    for (_iteratorAcrossCells += stride; _iteratorAcrossCells < _vectorOfCells->end(); _iteratorAcrossCells += stride) {
      if (not _iteratorAcrossCells->isEmpty() and isCellTypeBehaviorCorrect()) {
        _iteratorWithinOneCell = _iteratorAcrossCells->begin();
        break;
      }
    }
  }

  /**
   * checks if a cell has the correct cell type according to the behavior
   * @return true iff the cell type is proper according to the behavior
   */
  [[nodiscard]] bool isCellTypeBehaviorCorrect() const {
    // IMPORTANT: `this->` is necessary here! Without it clang 7, 8 and 9 fail due to an compiler bug:
    // https://stackoverflow.com/questions/55359614/clang-complains-about-constexpr-function-in-case-for-switch-statement
    switch (this->_behavior & ~IteratorBehavior::forceSequential) {
      case IteratorBehavior::ownedOrHaloOrDummy:
        [[fallthrough]];
      case IteratorBehavior::ownedOrHalo:
        return true;
      case IteratorBehavior::halo:
        return _flagManager->cellCanContainHaloParticles(_iteratorAcrossCells - _vectorOfCells->begin());
      case IteratorBehavior::owned:
        return _flagManager->cellCanContainOwnedParticles(_iteratorAcrossCells - _vectorOfCells->begin());
      default:
        utils::ExceptionHandler::exception(
            "ParticleIterator::isCellTypeBehaviorCorrect() encountered unknown iterator behavior: {}", this->_behavior);
        return false;
    }
  }

  /**
   * Indicates whether the particle has the correct owned state.
   * @return
   */
  [[nodiscard]] bool particleHasCorrectOwnershipState() const {
    // IMPORTANT: `this->` is necessary here! Without it clang 7, 8 and 9 fail due to an compiler bug:
    // https://stackoverflow.com/questions/55359614/clang-complains-about-constexpr-function-in-case-for-switch-statement
    if (this->_behavior == IteratorBehavior::ownedOrHaloOrDummy) {
      return true;
    } else if (_additionalParticleVectorToIterateState == AdditionalParticleVectorToIterateState::iterating) {
      return static_cast<unsigned int>(
                 (*_additionalVectors[_additionalVectorIndex])[_additionalVectorPosition].getOwnershipState()) &
             static_cast<unsigned int>(this->_behavior);
    } else {
      return static_cast<unsigned int>(_iteratorWithinOneCell->getOwnershipState()) &
             static_cast<unsigned int>(this->_behavior);
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
  IteratorBehavior _behavior;

  /**
   * Enum class that specifies the iterating state of the additional particle vector.
   */
  enum class AdditionalParticleVectorToIterateState { ignore, notStarted, iterating };

  /**
   * Specifies if this iterator should also iterate over the additional array.
   */
  AdditionalParticleVectorToIterateState _additionalParticleVectorToIterateState;

  /**
   * Vector of pointers to additional Particle Vectors this ParticleIterator will iterate over.
   */
  ParticleVecTypeInternal _additionalVectors;

  /**
   * Index of the additional vector that is currently processed.
   */
  size_t _additionalVectorIndex;

  /**
   * Index within the own additional particle vector.
   */
  size_t _additionalVectorPosition;
};
}  // namespace autopas::internal
