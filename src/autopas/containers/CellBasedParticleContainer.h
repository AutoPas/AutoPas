/**
 * @file CellBasedParticleContainer.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */
#pragma once

#include <algorithm>
#include <array>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

// consider multiple inheritance or delegation to avoid virtual call to Functor
/**
 * The CellBasedParticleContainer class stores particles in some object and provides
 * methods to iterate over its particles.
 * @tparam ParticleCell_T Class for the particle cells
 */
template <class ParticleCell_T>
class CellBasedParticleContainer : public ParticleContainerInterface<typename ParticleCell_T::ParticleType> {
  using ParticleType = typename ParticleCell_T::ParticleType;
  using ParticleCellType = ParticleCell_T;

 public:
  /**
   * Constructor of CellBasedParticleContainer
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   * @param sortingThreshold
   */
  CellBasedParticleContainer(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                             const double cutoff, const double skin, const size_t sortingThreshold)
      : ParticleContainerInterface<ParticleType>(skin),
        _cells(),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _cutoff(cutoff),
        _skin(skin),
        _sortingThreshold(sortingThreshold) {}

  /**
   * Destructor of CellBasedParticleContainer.
   */
  ~CellBasedParticleContainer() override = default;

  /**
   * Delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  CellBasedParticleContainer(const CellBasedParticleContainer &obj) = delete;

  /**
   * Delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  CellBasedParticleContainer &operator=(const CellBasedParticleContainer &other) = delete;

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const final { return _boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const final { return _boxMin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const final { return _cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) final { _cutoff = cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const final { return _cutoff + _skin; }
  /**
   * Returns the verlet Skin length
   * @return _skin
   */
  [[nodiscard]] double getVerletSkin() const final { return _skin; }

  /**
   * Deletes all particles from the container.
   */
  void deleteAllParticles() override {
    /// @todo: find a sensible value for magic number
    /// numThreads should be at least 1 and maximal max_threads
    AUTOPAS_OPENMP(parallel for num_threads(std::clamp(static_cast<int>(this->_cells.size()) / 1000,  \
                                                       1,                                             \
                                                       autopas::autopas_get_preferred_num_threads())))
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      this->_cells[i].clear();
    }
  }

  /**
   * @copydoc autopas::ParticleContainerInterface::getNumberOfParticles()
   */
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior) const override {
    size_t numParticles = 0ul;
    // parallelizing this loop is only worth it if we have LOTS of cells.
    // numThreads should be at least 1 and maximal max_threads
    AUTOPAS_OPENMP(parallel for num_threads(std::clamp(static_cast<int>(this->_cells.size()) / 100000, \
                                                       1,                                              \
                                                       autopas::autopas_get_preferred_num_threads()))  \
                                reduction(+ : numParticles))
    for (size_t index = 0; index < _cells.size(); ++index) {
      numParticles += _cells[index].getNumberOfParticles(behavior);
    }
    return numParticles;
  }

  /**
   * Get the total number of particles saved in the container (owned + halo + dummy).
   * @return Number of particles saved in the container (owned + halo + dummy).
   */
  [[nodiscard]] size_t size() const override {
    size_t numParticles = 0ul;
    // parallelizing this loop is only worth it if we have LOTS of cells.
    // numThreads should be at least 1 and maximal max_threads
    AUTOPAS_OPENMP(parallel for num_threads(std::clamp(static_cast<int>(this->_cells.size()) / 100000, \
                                                       1,                                              \
                                                       autopas::autopas_get_preferred_num_threads()))  \
                                reduction(+ : numParticles))
    for (size_t index = 0; index < _cells.size(); ++index) {
      numParticles += _cells[index].size();
    }
    return numParticles;
  }

  /**
   * Get immutable vector of cells.
   * @return immutable reference to _cells
   */
  [[nodiscard]] const std::vector<ParticleCellType> &getCells() const { return _cells; }

 protected:
  /**
   * Vector of particle cells.
   * All particle containers store their particles in ParticleCells. This is the
   * common vector for this purpose.
   */
  std::vector<ParticleCellType> _cells;
  /**
   * If the number of particles in a cell or cell pair exceeds this threshold, the particles will be sorted.
   * To be forwarded to cell traversals.
   */
  size_t _sortingThreshold;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
  double _skin;
};

}  // namespace autopas
