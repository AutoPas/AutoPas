/**
 * @file ParticleIterator.h
 *
 * @authors F. Gratl
 * @date 15.12.2022
 */

#pragma once

#include <array>
#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

// forward declaration
template <class Particle>
class ParticleContainerInterface;

namespace internal {
/**
 * Function to access private iterator.deleteCurrentParticle() via friend.
 * The function is private because it should not be called from outside of AutoPas since the logic handler has to keep
 * track of the number of particles.
 * This function is explicitly moved to the namespace "internal" so if you call it be sure you know what you do.
 *
 * @tparam ParticleIterator
 * @param iterator
 */
template <class ParticleIterator>
void deleteParticle(ParticleIterator &iterator) {
  iterator.deleteCurrentParticle();
}
}  // namespace internal

/**
 * Public iterator class that iterates over a particle container and additional vectors (which are typically stored in
 * the logic handler).
 * It supports parallelism over cells by being instantiated in a parallel region and particle deletion while iterating.
 * Inserting particles might invalidate the iterator.
 * There are no guarantees about the order in which particles are iterated.
 *
 * @tparam Particle
 * @tparam modifiable
 */
template <class Particle, bool modifiable>  // TODO: add template parameter for region iterator
class ContainerIterator {
  template <class T>
  friend void internal::deleteParticle(T &);

 public:
  /**
   * Type of the particle this iterator points to. Switch for const iterators.
   */
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;
  /**
   * Type of the additional vector collection. Switch for const iterators.
   */
  using ParticleVecType =
      std::conditional_t<modifiable, std::vector<std::vector<Particle> *>, std::vector<std::vector<Particle> const *>>;
  /**
   * Type of the Particle Container type. Switch for const iterators.
   */
  using ContainerType =
      std::conditional_t<modifiable, ParticleContainerInterface<Particle>, const ParticleContainerInterface<Particle>>;

  /**
   * Constructor meant to be called from the logic handler.
   * @note ATTENTION: This Iterator might be invalid after construction if there are no particles that satisfy it!
   *
   * @param container Reference to the particle container to iterate.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated over.
   * @param additionalVectorsToIterate Thread buffers of additional Particle vector to iterate over.
   */
  ContainerIterator(ContainerType &container, IteratorBehavior behavior, ParticleVecType *additionalVectorsToIterate)
      : _container(container),
        _behavior(behavior),
        _nextVectorIndex((behavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num()),
        _vectorIndexOffset((behavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads()) {
    if (additionalVectorsToIterate) {
      // store pointers to all additional vectors
      _additionalVectors.insert(_additionalVectors.end(), additionalVectorsToIterate->begin(),
                                additionalVectorsToIterate->end());
    }
    // fetches the next (=first) valid particle or sets _currentParticle = nullptr which marks the iterator as invalid.
    this->operator++();
  }

  /**
   * Increments the iterator.
   *
   * The idea is that this operator either queries the container with the "next" indices it got previously (or 0,0
   * initially), or if there is nothing left in the container it handles the iteration through the additional vectors.
   *
   * @return *this
   */
  inline ContainerIterator<Particle, modifiable> &operator++() {
    if (not _iteratingAdditionalVectors) {
      // this is not a region iter hence we stretch the bounding box to the numeric max
      constexpr std::array<double, 3> minBox{std::numeric_limits<double>::lowest(),
                                             std::numeric_limits<double>::lowest(),
                                             std::numeric_limits<double>::lowest()};
      constexpr std::array<double, 3> maxBox{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                                             std::numeric_limits<double>::max()};
      // this either gives us a particle with the desired properties or a nullptr
      std::tie(_currentParticle, _nextVectorIndex, _nextParticleIndex) =
          _container.getParticle(_nextVectorIndex, _nextParticleIndex, _behavior, minBox, maxBox);
      // if getParticle told us that the container doesn't have a particle in the first vector for our thread...
      if (_currentParticle == nullptr and not _additionalVectors.empty()) {
        // determine which additional vector this thread should start to iterate
        // we invert the assignment of threads to additional vectors because if there are more threads than cells
        // these surplus threads can then directly start with the additional vectors
        _nextVectorIndex = (_behavior & IteratorBehavior::forceSequential)
                               ? 0
                               : _additionalVectors.size() - 1 - autopas_get_thread_num();
        _nextParticleIndex = 0;
        _iteratingAdditionalVectors = true;
      } else {
        return *this;
      }
    }
    // no "else" here because the first case might trigger the second if we run to the end of the container
    if (_iteratingAdditionalVectors) {
      // check all vectors ...
      while (_nextVectorIndex < _additionalVectors.size()) {
        // ... and all particles within the vector ...
        while (_nextParticleIndex < _additionalVectors[_nextVectorIndex]->size()) {
          _currentParticle = &(_additionalVectors[_nextVectorIndex]->operator[](_nextParticleIndex));
          ++_nextParticleIndex;
          // ... until we find a particle that satisfies our requirements.
          if (particleHasCorrectOwnershipState(*_currentParticle)) {
            // if the particle is at the end of the current vector bump _nextVectorIndex and reset _nextParticleIndex.
            if (_nextParticleIndex >= _additionalVectors[_nextVectorIndex]->size()) {
              _nextVectorIndex += _vectorIndexOffset;
              _nextParticleIndex = 0;
            }
            return *this;
          }
        }
        _nextVectorIndex += _vectorIndexOffset;
        _nextParticleIndex = 0;
      }
    }
    // if we reach this point there is no satisfying particle left neither in the container nor the additional vectors.
    _currentParticle = nullptr;
    return *this;
  }

  /**
   * Dereference operator.
   * @return Reference to the current particle.
   */
  inline ParticleType &operator*() const { return *_currentParticle; }

  /**
   * Dereference operator.
   * @return Pointer to the current particle.
   */
  inline ParticleType *operator->() const { return _currentParticle; }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  [[nodiscard]] bool isValid() const { return _currentParticle != nullptr; }

  /**
   * Checks if the current iterator has a given validity.
   * @param input
   * @return
   */
  bool operator==(const bool input) const { return isValid() == input; }

  /**
   * Checks if the current iterator does not have a given validity.
   * @param input
   * @return
   */
  bool operator!=(const bool input) const { return not(*this == input); }

 private:
  /**
   * Indicates whether the particle has the correct owned state.
   * @return True iff the given particle's ownership state is compatible with the iterator behavior.
   */
  [[nodiscard]] bool particleHasCorrectOwnershipState(const Particle &p) const {
    // Dummy is not in sync with iterator behavior because it needs to be 0.
    if (_behavior & IteratorBehavior::ownedOrHaloOrDummy or (_behavior == IteratorBehavior::dummy and p.isDummy())) {
      return true;
    } else {
      // relies on both enums having the same encoding. This should be guaranteed statically in IteratorBehaviorTest!
      return static_cast<unsigned int>(p.getOwnershipState()) & static_cast<unsigned int>(_behavior);
    }
  }

  /**
   * Deletes the particle the particle currently pointed to. This function uses swap-delete and thus will change
   * the order of elements in the container. After deletion _currentParticle points to the next particle,
   * which is the particle that was swapped here. So when used in a loop do not increment the iterator after deletion.
   */
  void deleteCurrentParticle() {
    if (_iteratingAdditionalVectors) {
      // the current particle address needs to be between start and end of the current vector
      // otherwise it is from the previous vector
      const auto currentVectorIndex = (_currentParticle >= &_additionalVectors[_nextVectorIndex]->front() and
                                       _currentParticle <= &_additionalVectors[_nextVectorIndex]->back())
                                          ? _nextVectorIndex
                                          : _nextVectorIndex - _vectorIndexOffset;
      auto &currentVector = *_additionalVectors[currentVectorIndex];
      // swap-delete
      *_currentParticle = currentVector.back();
      currentVector.pop_back();
      // make sure _currentParticle always points to a valid particle if there is one.
      // FIXME: region iter
      if (currentVector.empty() or not particleHasCorrectOwnershipState(*_currentParticle)) {
        this->operator++();
      }
    } else {
      const auto pointerValid = _container.deleteParticle(*_currentParticle);
      // FIXME: region iter
      if (not pointerValid or not particleHasCorrectOwnershipState(*_currentParticle)) {
        this->operator++();
      }
    }
  }

  /**
   * Pointer to container that is iterated.
   */
  ContainerType &_container;

  /**
   * Indicator which type of particles this iterator covers.
   */
  IteratorBehavior _behavior;

  /**
   * Flag indicating if we are currently iterating the additional vectors.
   * Since we always check the container first, this starts as false.
   */
  bool _iteratingAdditionalVectors{false};

  /**
   * Vector of pointers to additional Particle Vectors this ParticleIterator will iterate over.
   */
  ParticleVecType _additionalVectors;

  /**
   * The particle the iterator currently points to.
   * This is always either a valid particle satisfying all iterator requirements or nullptr.
   */
  ParticleType *_currentParticle = nullptr;

  /**
   * Index of the Vector where the next particle is found. This might be the same index as the current one.
   * "Vector" typically refers to either a cell in the particle container or one of the additional vectors.
   */
  size_t _nextVectorIndex;

  /**
   * Offset which which to iterate over vectors. Determined through number of threads used for iterating.
   */
  size_t _vectorIndexOffset;

  /**
   * Index within the next vector
   */
  size_t _nextParticleIndex{0};
};
}  // namespace autopas
