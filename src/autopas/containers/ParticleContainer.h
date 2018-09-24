/**
 * @file ParticleContainer.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversalInterface.h"
#include "autopas/pairwiseFunctors/Functor.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

namespace autopas {

// consider multiple inheritance or delegation to avoid virtual call to Functor
/**
 * The ParticleContainer class stores particles in some object and provides
 * methods to iterate over its particles.
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle, class ParticleCell, class SoAArraysType = typename Particle::SoAArraysType>
class ParticleContainer : public ParticleContainerInterface<Particle, ParticleCell> {
 private:
  static const std::vector<TraversalOptions> &DefaultApplicableTraversals() {
    static const std::vector<TraversalOptions> v{};
    return v;
  }

 public:
  /// type of the Particle
  typedef Particle ParticleType;

  /// type of the ParticleCell
  typedef ParticleCell ParticleCellType;
  /**
   * Constructor of ParticleContainer
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param applicableTraversals Traversals applicable for this Container
   */
  ParticleContainer(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                    const std::vector<TraversalOptions> &applicableTraversals = DefaultApplicableTraversals())
      : _cells(), _applicableTraversals(applicableTraversals), _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff) {}

  /**
   * destructor of ParticleContainer
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
   * Get the upper corner of the container
   * @return upper corner of the container
   */
  const std::array<double, 3> &getBoxMax() const override final { return _boxMax; }

  /**
   * Set the upper corner of the container
   * @param boxMax upper corner to be set
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final { _boxMax = boxMax; }

  /**
   * Get the lower corner of the container
   * @return lower corner of the container
   */
  const std::array<double, 3> &getBoxMin() const override final { return _boxMin; }

  /**
   * Set the lower corner of the container
   * @param boxMin lower corner to be set
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final { _boxMin = boxMin; }

  /**
   * Return the cutoff of the container
   * @return
   */
  double getCutoff() const override final { return _cutoff; }

  /**
   * Set the cutoff of the container
   * @param cutoff
   */
  void setCutoff(double cutoff) override final { _cutoff = cutoff; }

  /**
   * Checks if the given traversals are applicable to this traversal.
   * @param traversalOptions
   * @return True iff traversalOptions is a subset of _applicableTraversals
   */
  bool checkIfTraversalsAreApplicable(std::vector<TraversalOptions> traversalOptions) {
    for (auto &option : traversalOptions) {
      if (find(_applicableTraversals.begin(), _applicableTraversals.end(), option) == _applicableTraversals.end())
        return false;
    }
    return true;
  }

  /**
   * Deletes all particles from the container.
   */
  void deleteAllParticles() override {
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for magic number
    // numThreads should be at least 1 and maximal max_threads
    int numThreads = std::max(1, std::min(omp_get_max_threads(), (int)(this->_cells.size() / 1000)));
    AutoPasLog(trace, "Using {} threads", numThreads);
#pragma omp parallel for num_threads(numThreads)
#endif
    for (size_t i = 0; i < this->_cells.size(); ++i) {
      this->_cells[i].clear();
    }
  }

  unsigned long getNumParticles() override {
    size_t numParticles = 0ul;
#ifdef AUTOPAS_OPENMP
    // @todo: find a sensible value for magic number
    // numThreads should be at least 1 and maximal max_threads
    int numThreads = std::max(1, std::min(omp_get_max_threads(), (int)(this->_cells.size() / 1000)));
    AutoPasLog(trace, "Using {} threads", numThreads);
#pragma omp parallel for num_threads(numThreads) reduction(+ : numParticles)
#endif
    for (size_t index = 0; index < _cells.size(); ++index) {
      numParticles += _cells[index].numParticles();
    }
    return numParticles;
  }

 protected:
  /**
   * Vector of particle cells.
   * All particle containers store their particles in ParticleCells. This is the
   * common vector for this purpose.
   */
  std::vector<ParticleCell> _cells;
  /**
   * Vector of all applicable traversal options for the container.
   */
  const std::vector<TraversalOptions> &_applicableTraversals;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
};

}  // namespace autopas
