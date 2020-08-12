//
// Created by lunaticcoding on 04.05.20.
//

#pragma AUTOPAS_PARTICLELIST_H

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
  using particleListImpType = std::vector<Type>;

  using iterator = typename particleListImpType::iterator;
  using const_iterator = typename particleListImpType::const_iterator;

 public:
  ParticleVector<Type>() {
    _dirty = false;
    _dirtyIndex = 0;
    particleListImp = std::vector<Type>();
  }

  /**
   * Returns the dirty flag, indicating whether some of the Particle references are out of date.
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
   * Add a Particle to the data structure.
   * @param value A reference to the value to be stored
   */
  void push_back(Type &value) {
    particleListLock.lock();
    _dirty |= particleListImp.capacity() == particleListImp.size();
    if (_dirty) {
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
   * Begin of the iterator over dirty Particles
   * @return Start of the iterator
   */
  iterator beginDirty() { return particleListImp.begin() + _dirtyIndex; }
  /**
   * End of the iterator over dirty Particles
   * @return End of the iterator
   */
  iterator endDirty() { return particleListImp.end(); }

 private:
  bool _dirty;
  int _dirtyIndex;
  autopas::AutoPasLock particleListLock;
  std::vector<Type> particleListImp;
};
