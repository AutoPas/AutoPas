/**
 * @file ParticleContainer.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <selectors/TraversalSelector.h>
#include <array>
#include "ParticleContainerInterface.h"
#include "pairwiseFunctors/Functor.h"

namespace autopas {

// consider multiple inheritance or delegation to avoid virtual call to Functor
/**
 * The ParticleContainer class stores particles in some object and provides
 * methods to iterate over its particles.
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle, class ParticleCell, class SoAArraysType = typename Particle::SoAArraysType>
class ParticleContainer : public ParticleContainerInterface<Particle> {
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
      : _cells(),
        _traversalSelector(nullptr),    // needs to be instantiated by respective container.
        _applicableTraversals(applicableTraversals),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _cutoff(cutoff) {}

  /**
   * destructor of ParticleContainer
   */
  ~ParticleContainer() override = default;

  /**
   * delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  ParticleContainer(const ParticleContainer &obj) = delete;

  /**
   * delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  ParticleContainer &operator=(const ParticleContainer &other) = delete;

  /**
   * get the upper corner of the container
   * @return upper corner of the container
   */
  const std::array<double, 3> &getBoxMax() const override final { return _boxMax; }

  /**
   * set the upper corner of the container
   * @param boxMax upper corner to be set
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final { _boxMax = boxMax; }

  /**
   * get the lower corner of the container
   * @return lower corner of the container
   */
  const std::array<double, 3> &getBoxMin() const override final { return _boxMin; }

  /**
   * set the lower corner of the container
   * @param boxMin lower corner to be set
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final { _boxMin = boxMin; }

  /**
   * return the cutoff of the container
   * @return
   */
  double getCutoff() const override final { return _cutoff; }

  /**
   * set the cutoff of the container
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
   * Determine the optimal traversal for the current situation.
   * @tparam PairwiseFunctor
   * @tparam useSoA
   * @tparam useNewton3
   * @param pairwiseFunctor Functor to optimize for.
   */
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  void tuneTraversal(PairwiseFunctor &pairwiseFunctor) {
    if (_traversalSelector != nullptr)
      _traversalSelector->template tune<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
  }

 protected:
  /**
   * vector of particle cells.
   * All particle containers store their particles in ParticleCells. This is the
   * common vector for this purpose.
   */
  std::vector<ParticleCell> _cells;
  /**
   * Selector for traversal of the container.
   */
  std::unique_ptr<TraversalSelector<ParticleCell>> _traversalSelector;
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
