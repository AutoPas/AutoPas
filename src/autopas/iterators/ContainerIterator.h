/**
 * @file ContainerIterator.h
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

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

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

namespace containerIteratorUtils {
/**
 * Indicates whether the particle has the correct ownership state and if this is a region iterator is in the region.
 *
 * @tparam regionIter
 * @tparam Particle
 * @tparam Arr Template because if this is a region iterator ContainerIterator will pass an empty struct
 * @param p
 * @param behavior
 * @param regionMin If regionIter is true this should be a std::array<double, 3>
 * @param regionMax If regionIter is true this should be a std::array<double, 3>
 * @return True iff the given particle's ownership state is compatible with the iterator behavior.
 */
template <bool regionIter, class Particle, class Arr>
[[nodiscard]] bool particleFulfillsIteratorRequirements(const Particle &p, IteratorBehavior behavior,
                                                        const Arr &regionMin, const Arr &regionMax) {
  bool particleOk = false;
  // Check ownership
  // Dummy is not in sync with iterator behavior because it needs to be 0.
  // `a & b == b` idiom is to check a == b disregarding any bits beyond b. This is needed due to e.g. forceSequential.
  if (((behavior & IteratorBehavior::ownedOrHaloOrDummy) == IteratorBehavior::ownedOrHaloOrDummy) or
      (behavior & IteratorBehavior::dummy and p.isDummy())) {
    particleOk = true;
  } else {
    // relies on both enums having the same encoding. This should be guaranteed statically in IteratorBehaviorTest!
    particleOk = static_cast<unsigned int>(p.getOwnershipState()) & static_cast<unsigned int>(behavior);
  }
  if constexpr (regionIter) {
    // If this is a region iterator check if the particle is in the box.
    if (particleOk and not utils::inBox(p.getR(), regionMin, regionMax)) {
      particleOk = false;
    }
  }
  return particleOk;
}
}  // namespace containerIteratorUtils

/**
 * Public iterator class that iterates over a particle container and additional vectors (which are typically stored in
 * the logic handler).
 * It supports parallelism over cells by being instantiated in a parallel region and particle deletion while iterating.
 * Inserting particles might invalidate the iterator.
 * There are no guarantees about the order in which particles are iterated.
 *
 * @tparam Particle
 * @tparam modifiable If false, this is a const iterator.
 * @tparam regionIter If false, avoid any region checks and iterate the whole container.
 */
template <class Particle, bool modifiable, bool regionIter>
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
   * Default constructor that guarantees an invalid iterator.
   */
  ContainerIterator() : _currentParticle(nullptr){};

  /**
   * Region Iterator constructor meant to be called from the logic handler.
   *
   * @param container Reference to the particle container to iterate.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated over.
   * @param additionalVectorsToIterate Thread buffers of additional Particle vector to iterate over.
   * @param regionMin Left Front Lower corner of the iterator's region.
   * @param regionMax Right Back Upper corner of the iterator's region.
   */
  ContainerIterator(ContainerType &container, IteratorBehavior behavior, ParticleVecType *additionalVectorsToIterate,
                    const std::array<double, 3> &regionMin, const std::array<double, 3> &regionMax)
      : ContainerIterator(nullptr, container, behavior, additionalVectorsToIterate, regionMin, regionMax) {
    // sanity check
    static_assert(regionIter == true,
                  "Constructor for Region iterator called but template argument regionIter is false");
  }

  /**
   * Regular Iterator constructor meant to be called from the logic handler.
   *
   * @param container Reference to the particle container to iterate.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated over.
   * @param additionalVectorsToIterate Thread buffers of additional Particle vector to iterate over.
   */
  ContainerIterator(ContainerType &container, IteratorBehavior behavior, ParticleVecType *additionalVectorsToIterate)
      : ContainerIterator(nullptr, container, behavior, additionalVectorsToIterate, {}, {}) {
    static_assert(regionIter == false,
                  "Constructor for non-Region iterator called but template argument regionIter is true");
  }

  /**
   * Copy constructor
   * @param other
   */
  ContainerIterator(const ContainerIterator<Particle, modifiable, regionIter> &other)
      : _container(other._container),
        _currentParticleIndex(other._currentParticleIndex),
        _currentVectorIndex(other._currentVectorIndex),
        _currentParticle(other._currentParticle),
        _additionalVectors(other._additionalVectors),
        _behavior(other._behavior),
        _iteratingAdditionalVectors(other._iteratingAdditionalVectors),
        _vectorIndexOffset(other._vectorIndexOffset),
        _regionMin(other._regionMin),
        _regionMax(other._regionMax) {}

  /**
   * Copy assignment operator
   * @param other
   * @return
   */
  ContainerIterator<Particle, modifiable, regionIter> &operator=(
      const ContainerIterator<Particle, modifiable, regionIter> &other) {
    if (this != other) {
      _container = other._container;
      _currentParticleIndex = other._currentParticleIndex;
      _currentVectorIndex = other._currentVectorIndex;
      _currentParticle = other._currentParticle;
      _additionalVectors = other._additionalVectors;
      _behavior = other._behavior;
      _iteratingAdditionalVectors = other._iteratingAdditionalVectors;
      _vectorIndexOffset = other._vectorIndexOffset;
      if constexpr (regionIter) {
        _regionMin = other._regionMin;
        _regionMax = other._regionMax;
      }
    }
    return *this;
  }

  /**
   * Move constructor
   * @param other
   */
  ContainerIterator(ContainerIterator<Particle, modifiable, regionIter> &&other) noexcept
      : _container(other._container),
        _currentParticleIndex(other._currentParticleIndex),
        _currentVectorIndex(other._currentVectorIndex),
        _currentParticle(other._currentParticle),
        _additionalVectors(std::move(other._additionalVectors)),
        _behavior(other._behavior),
        _iteratingAdditionalVectors(other._iteratingAdditionalVectors),
        _vectorIndexOffset(other._vectorIndexOffset),
        _regionMin(std::move(other._regionMin)),
        _regionMax(std::move(other._regionMax)) {}

  /**
   * Move assignment operator
   * @param other
   * @return
   */
  ContainerIterator<Particle, modifiable, regionIter> &operator=(
      ContainerIterator<Particle, modifiable, regionIter> &&other) noexcept {
    if (this != &other) {
      _container = other._container;
      _currentParticleIndex = other._currentParticleIndex;
      _currentVectorIndex = other._currentVectorIndex;
      _currentParticle = other._currentParticle;
      _additionalVectors = std::move(other._additionalVectors);
      _behavior = other._behavior;
      _iteratingAdditionalVectors = other._iteratingAdditionalVectors;
      _vectorIndexOffset = other._vectorIndexOffset;
      if constexpr (regionIter) {
        _regionMin = std::move(other._regionMin);
        _regionMax = std::move(other._regionMax);
      }
    }
    return *this;
  }

 private:
  /**
   * Actual internal constructor.
   *
   * The dummy parameter exists to distinguish this constructor's signature from the others.
   *
   * @param container Reference to the particle container to iterate.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated over.
   * @param additionalVectorsToIterate Thread buffers of additional Particle vector to iterate over.
   */
  ContainerIterator(void * /*dummy*/, ContainerType &container, IteratorBehavior behavior,
                    ParticleVecType *additionalVectorsToIterate, const std::array<double, 3> &regionMin,
                    const std::array<double, 3> &regionMax)
      : _container(&container),
        _behavior(behavior),
        _currentVectorIndex(0),
        _vectorIndexOffset((behavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads()) {
    if (additionalVectorsToIterate) {
      // store pointers to all additional vectors
      _additionalVectors.insert(_additionalVectors.end(), additionalVectorsToIterate->begin(),
                                additionalVectorsToIterate->end());
    }

    if constexpr (regionIter) {
      // clamp region to the container, either with or without halo
      const auto boxMax = (_behavior & IteratorBehavior::halo)
                              ? utils::ArrayMath::addScalar(container.getBoxMax(), container.getInteractionLength())
                              : container.getBoxMax();
      const auto boxMin = (_behavior & IteratorBehavior::halo)
                              ? utils::ArrayMath::subScalar(container.getBoxMin(), container.getInteractionLength())
                              : container.getBoxMin();
      _regionMax = utils::ArrayMath::min(regionMax, boxMax);
      _regionMin = utils::ArrayMath::max(regionMin, boxMin);
    }
    // fetches the next (=first) valid particle or sets _currentParticle = nullptr which marks the iterator as invalid.
    fetchParticleAtCurrentIndex();
  }

 public:
  /**
   * Increments the iterator.
   *
   * The idea is that this operator either queries the container with the current indices, or, if there is nothing left
   * in the container, it handles the iteration through the additional vectors.
   *
   * @return *this
   */
  inline ContainerIterator<Particle, modifiable, regionIter> &operator++() {
    // bump the index. If it is invalid now, the container will give us the proper one.
    ++_currentParticleIndex;
    fetchParticleAtCurrentIndex();
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
   * This function queries the container with the current iterator indices and updates them as well as the particle
   * pointer according to the container's response.
   */
  void fetchParticleAtCurrentIndex() {
    if (not _iteratingAdditionalVectors) {
      // getParticle either gives us a particle with the desired properties or a nullptr
      if constexpr (regionIter) {
        std::tie(_currentParticle, _currentVectorIndex, _currentParticleIndex) =
            _container->getParticle(_currentVectorIndex, _currentParticleIndex, _behavior, _regionMin, _regionMax);
      } else {
        std::tie(_currentParticle, _currentVectorIndex, _currentParticleIndex) =
            _container->getParticle(_currentVectorIndex, _currentParticleIndex, _behavior);
      }
      // if getParticle told us that the container doesn't have a particle in the first vector for our thread...
      if (_currentParticle == nullptr and not _additionalVectors.empty()) {
        // determine which additional vector this thread should start to iterate
        // we invert the assignment of threads to additional vectors because if there are more threads than cells
        // these surplus threads can then directly start with the additional vectors
        _currentVectorIndex = (_behavior & IteratorBehavior::forceSequential)
                                  ? 0
                                  : autopas_get_num_threads() - 1 - autopas_get_thread_num();
        _currentParticleIndex = 0;
        _iteratingAdditionalVectors = true;
      } else {
        // Case: nothing left in the container and no additional vectors to iterate
        return;
      }
    }
    // no "else" here because the first case might trigger the second if we reached the end of the container
    if (_iteratingAdditionalVectors) {
      // check all vectors ...
      for (; _currentVectorIndex < _additionalVectors.size(); _currentVectorIndex += _vectorIndexOffset) {
        // ... and all particles within the vector ...
        for (; _currentParticleIndex < _additionalVectors[_currentVectorIndex]->size(); ++_currentParticleIndex) {
          _currentParticle = &(_additionalVectors[_currentVectorIndex]->operator[](_currentParticleIndex));
          // ... until we find a particle that satisfies our requirements.
          if (containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(*_currentParticle, _behavior,
                                                                                       _regionMin, _regionMax)) {
            return;
          }
        }
        _currentParticleIndex = 0;
      }
    }
    // if we reach this point there is no satisfying particle left neither in the container nor the additional vectors.
    _currentParticle = nullptr;
    _currentParticleIndex = std::numeric_limits<decltype(_currentParticleIndex)>::max();
    _currentVectorIndex = std::numeric_limits<decltype(_currentVectorIndex)>::max();
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
      auto &currentVector = *_additionalVectors[_currentVectorIndex];
      // swap-delete
      *_currentParticle = currentVector.back();
      currentVector.pop_back();
      // make sure _currentParticle always points to a valid particle if there is one.
      if (currentVector.empty() or not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
                                       *_currentParticle, _behavior, _regionMin, _regionMax)) {
        this->operator++();
      }
    } else {
      const auto indicesValid = _container->deleteParticle(_currentVectorIndex, _currentParticleIndex);
      // Edge cases:
      if (not indicesValid) {
        // CASE: the indices and thus the pointer are invalid now
        this->operator++();
      } else if (not containerIteratorUtils::particleFulfillsIteratorRequirements<regionIter>(
                     *_currentParticle, _behavior, _regionMin, _regionMax)) {
        // CASE: the particle was swapped and now the pointer points at something uninteresting
        this->operator++();
      }
    }
  }

  /**
   * Pointer to container that is iterated.
   */
  ContainerType *_container;

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
   * Index of the Vector where the current particle is found.
   * "Vector" typically refers to either a cell in the particle container or one of the additional vectors.
   */
  size_t _currentVectorIndex{};

  /**
   * Offset which which to iterate over vectors. Determined through number of threads used for iterating.
   */
  size_t _vectorIndexOffset{};

  /**
   * Index within the current vector.
   */
  size_t _currentParticleIndex{0};

  /**
   * Dummy type for a data type of size zero.
   */
  struct empty {};
  /**
   * Depending on whether this is a region iterator define the type of the members to be three doubles or nothing.
   */
  using RegionCornerT = std::conditional_t<regionIter, std::array<double, 3>, empty>;
  /**
   * Lower corner of the region iterator box.
   * @note Marked as [[no_unique_address]] so that regular iterator objects don't allocate memory for these unnecessary
   * fields.
   */
  [[no_unique_address]] RegionCornerT _regionMin{};
  /**
   * Upper corner of the region iterator box.
   * @note Marked as [[no_unique_address]] so that regular iterator objects don't allocate memory for these unnecessary
   * fields.
   */
  [[no_unique_address]] RegionCornerT _regionMax{};
};
}  // namespace autopas
