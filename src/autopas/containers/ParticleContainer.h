/**
 * @file ParticleContainer.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>

#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/TraversalInterface.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

namespace autopas {

// consider multiple inheritance or delegation to avoid virtual call to Functor
/**
 * The ParticleContainer class stores particles in some object and provides
 * methods to iterate over its particles.
 * @tparam ParticleCell Class for the particle cells
 */
template <class ParticleCell, class SoAArraysType = typename ParticleCell::ParticleType::SoAArraysType>
class ParticleContainer : public ParticleContainerInterface<typename ParticleCell::ParticleType> {
 public:
  /**
   * Constructor of ParticleContainer
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param skin
   */
  ParticleContainer(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                    const double skin)
      : _cells(), _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff), _skin(skin) {}

      using ParticleCellType = ParticleCell;

    /**
     * Destructor of ParticleContainer.
     */
    ParticleCellTypeEnum getParticleCellTypeEnum() override  {
        ParticleCell someCell = ParticleCell();
        return someCell.getParticleCellTypeAsEnum();
    };

  /**
   * Destructor of ParticleContainer.
   */
  ~ParticleContainer() override = default;

  /**
   * Delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  ParticleContainer(const ParticleContainer &obj) = delete;

  /**
   * Delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  ParticleContainer &operator=(const ParticleContainer &other) = delete;

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMax()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const override final { return _boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMax()
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final { _boxMax = boxMax; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getBoxMin()
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const override final { return _boxMin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setBoxMin()
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final { _boxMin = boxMin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getCutoff()
   */
  [[nodiscard]] double getCutoff() const override final { return _cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setCutoff()
   */
  void setCutoff(double cutoff) override final { _cutoff = cutoff; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getSkin()
   */
  [[nodiscard]] double getSkin() const override final { return _skin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::setSkin()
   */
  void setSkin(double skin) override final { _skin = skin; }

  /**
   * @copydoc autopas::ParticleContainerInterface::getInteractionLength()
   */
  [[nodiscard]] double getInteractionLength() const override final { return _cutoff + _skin; }

  /**
   * Deletes all particles from the container.
   */
  void deleteAllParticles() override {
#ifdef AUTOPAS_OPENMP
    /// @todo: find a sensible value for magic number
    /// numThreads should be at least 1 and maximal max_threads
    int numThreads = std::max(1, std::min(omp_get_max_threads(), (int)(this->_cells.size() / 1000)));
    AutoPasLog(trace, "Using {} threads", numThreads);
#pragma omp parallel for num_threads(numThreads)
#endif
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      this->_cells[i].clear();
    }
  }

  /**
   * Get the number of particles saved in the container.
   * @return Number of particles in the container.
   */
  [[nodiscard]] unsigned long getNumParticles() const override {
    size_t numParticles = 0ul;
#ifdef AUTOPAS_OPENMP
    /// @todo: find a sensible value for magic number
    /// numThreads should be at least 1 and maximal max_threads
    int numThreads = std::max(1, std::min(omp_get_max_threads(), (int)(this->_cells.size() / 1000)));
    AutoPasLog(trace, "Using {} threads", numThreads);
#pragma omp parallel for num_threads(numThreads) reduction(+ : numParticles)
#endif
    for (size_t index = 0; index < _cells.size(); ++index) {
      numParticles += _cells[index].numParticles();
    }
    return numParticles;
  }

  /**
   * Get immutable vector of cells.
   * @return immutable reference to _cells
   */
  const std::vector<ParticleCell> &getCells() const { return _cells; }

 protected:
  /**
   * Vector of particle cells.
   * All particle containers store their particles in ParticleCells. This is the
   * common vector for this purpose.
   */
  std::vector<ParticleCell> _cells;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
  double _skin;
};

}  // namespace autopas
