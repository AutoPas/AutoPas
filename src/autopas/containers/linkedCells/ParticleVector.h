/**
 * @file ParticleVector.h
 *
 * @author lunaticcoding
 * @date 04.05.2020
 */

#pragma once

#include <vector>

#include "autopas/utils/WrapOpenMP.h"

/**
 * ParticleVector class.
 * This class uses a std::vector to store particles and allows access to its references.
 * It also keeps track, which references are no longer valid, which happens
 * when operations like delete and resize are performed.
 * @tparam Type of the Particle that is being stored
 */
template <class Type>
class ParticleVector {
 public:
  ParticleVector<Type>() = default;

  /**
   * Returns the dirty flag, indicating whether Particles exist in the vector that are not stored in a cell yet.
   * @return True if dirty, false otherwise
   */
  bool isDirty() { return _dirty; }

  /**
   * Marks the ParticleVector as clean. Should be called after updating the references of dirty particles.
   */
  void markAsClean() {
    _dirty = false;
    _dirtyIndex = particleListImp.size();
  }

  /**
   * Remove all halo particles from the container and mark it as dirty.
   */
  void clearHaloParticles() {
    particleListImp.erase(std::remove_if(particleListImp.begin(), particleListImp.end(),
                                         [](const auto &particle) { return particle.isHalo(); }),
                          particleListImp.end());
    _dirty = true;
    _dirtyIndex = 0;
  }

  /**
   * Remove all dummy particles from the container and mark it as dirty.
   */
  void deleteDummyParticles() {
    particleListImp.erase(std::remove_if(particleListImp.begin(), particleListImp.end(),
                                         [](const auto &particle) { return particle.isDummy(); }),
                          particleListImp.end());
    _dirty = true;
    _dirtyIndex = 0;
  }

  /**
   * Add a Particle to the data structure.
   * @param value A reference to the value to be stored
   */
  void push_back(Type &value) {
    particleListLock.lock();
    _dirty = true;
    if (particleListImp.capacity() == particleListImp.size()) {
      _dirtyIndex = 0;
    }
    particleListImp.push_back(value);
    particleListLock.unlock();
  }

  /**
   * Get the number of Particles in the data structure.
   * @return Total number of Particles
   */
  int totalSize() { return particleListImp.size(); }

  /**
   * Get the number of dirty Particles in the data structure.
   * @return Number of dirty Particles
   */
  int dirtySize() { return totalSize() - _dirtyIndex; }

  /**
   * Indicates, whether References already stored in cells need to be updated.
   * @return Boolean indicating the above
   */
  bool needsRebuild() { return _dirtyIndex == 0; }

  /**
   * Begin of the iterator over dirty Particles
   * @return Start of the iterator
   */
  auto beginDirty() { return particleListImp.begin() + _dirtyIndex; }
  /**
   * End of the iterator over dirty Particles
   * @return End of the iterator
   */
  auto endDirty() { return particleListImp.end(); }

 private:
  /**
   * Flag indicating whether there are out-of-date references in the vector.
   */
  bool _dirty = false;
  /**
   * Index of the first out-of-date reference.
   */
  size_t _dirtyIndex{0};
  autopas::AutoPasLock particleListLock;
  std::vector<Type> particleListImp;
};
