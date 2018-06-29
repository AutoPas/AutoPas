/**
 * @file SingleCellIteratorInterface.h
 *
 * @date 31 May 2018
 * @author seckler
 */

#pragma once
#include "ParticleIteratorInterface.h"

namespace autopas {

/**
 * SingleCellIteratorInterface class to loop over particles of a single cell.
 *
 * @tparam Particle type of the Particles
 */
template <class Particle>
class SingleCellIteratorInterface : public ParticleIteratorInterface<Particle> {
 public:
  /**
   * default constructor of SingleCellIteratorInterface
   */
  SingleCellIteratorInterface() = default;

  /**
   * destructor of SingleCellIteratorInterface
   */
  virtual ~SingleCellIteratorInterface() = default;

  /**
   * equality operator.
   * if both iterators are invalid or if they point to the same particle, this
   * returns true
   * @param rhs
   * @return
   */
  virtual bool operator==(const SingleCellIteratorInterface<Particle> &rhs) const = 0;

  /**
   * inequality operator
   * descrition see operator==
   * @param rhs
   * @return
   */
  virtual bool operator!=(const SingleCellIteratorInterface<Particle> &rhs) const = 0;

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  virtual size_t getIndex() const = 0;
};

namespace internal {

/**
 * All implementations of the interface should inherit from this class. It extends the interface just by the clone
 * method, which is needed by the Wrapper.
 * @tparam Particle
 */
template <class Particle>
class SingleCellIteratorInterfaceImpl : public SingleCellIteratorInterface<Particle> {
 public:
  /**
   * Clones the current object, should allocate new object and return it.
   * @return the clone
   */
  virtual SingleCellIteratorInterfaceImpl<Particle> *clone() const = 0;
};
}  // namespace internal
}  // namespace autopas

/**
 * gets a static cell iterator from an iteratorwrapper of a cell.
 * @param iter the iterator to be defined
 * @param cell the cell for which the iterator should be get
 * @param body the body to be executed with the static iterator
 */
#define AUTOPAS_WITH_STATIC_CELL_ITER(iter, cell, body)                                                               \
  auto __wrapper = cell.begin();                                                                                      \
  auto __ptr = __wrapper.get();                                                                                       \
  {                                                                                                                   \
    if (auto __##iter##ptr = dynamic_cast<                                                                            \
            internal::SingleCellIterator<Particle, typename std::remove_reference<decltype(cell)>::type> *>(__ptr)) { \
      auto iter = *__##iter##ptr;                                                                                     \
      body                                                                                                            \
    } else if (auto __##iter##ptr = dynamic_cast<RMMParticleCellIterator<Particle> *>(__ptr)) {                       \
      auto iter = *__##iter##ptr;                                                                                     \
      body                                                                                                            \
    } else {                                                                                                          \
      autopas::utils::ExceptionHandler::exception("unknown iteratortype in WITH_STATIC_CELL_ITER");                   \
    }                                                                                                                 \
  }
