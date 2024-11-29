/**
 * @file HierarchicalGrid.h
 *
 * @date 29 Nov 2024
 * @author atacann
 */
#pragma once

#include <algorithm>
#include <array>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBlock3D.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/containers/cellTraversals/BalancedTraversal.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/iterators/ContainerIterator.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

/**
 * HierarchicalGrid class.
 * This class stores multiple LinkedCells containers, each for different sizes of particles
 * Traverses them all one by one, and cross hierarchy afterward
 * @tparam Particle type of the Particle
 */
template <class Particle>
class HierarchicalGrid : public ParticleContainerInterface<Particle> {
 public:
  /**
   *  Type of the ParticleCell.
   */
  using ParticleCell = FullParticleCell<Particle>;

  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle;

    /**
   * Constructor of the HierarchicalGrid class
   * @param boxMin
   * @param boxMax
   * @param baseCutoff base cutoff for each particle, it will be scaled by the size of a Particle
   * @param cutoffs cutoffs for each level of the hierarchy
   * @param skinPerTimestep
   * @param rebuildFrequency
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator the load estimation algorithm for balanced traversals.
   */
  HierarchicalGrid(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax, const double baseCutoff, 
              std::vector<double> cutoffs, const double skinPerTimestep, const unsigned int rebuildFrequency, const double cellSizeFactor = 1.0,
              LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell)
      : ParticleContainerInterface<Particle>(skinPerTimestep),
        _boxMin(boxMin),
        _boxMax(boxMax),
        _baseCutoff(baseCutoff),
        _cutoffs(cutoffs),
        _skin(skinPerTimestep * rebuildFrequency),
        _numHierarchyLevels(cutoffs.size()),
        _loadEstimator(loadEstimator) {
          if(_cutoffs.size() == 0) {
            AutoPasLog(ERROR, "Hierarchical Grid cutoffs vector is empty.");
            utils::ExceptionHandler::exception("Error: Hierarchical Grid cutoffs vector is empty.");
          }
          // make sure cutoffs are sorted
          std::sort(_cutoffs.begin(), _cutoffs.end());
          // generate LinkedCells for each hierarchy, with different cutoffs
          _hierarchies.reserve(_numHierarchyLevels);
          for(size_t i = 0; i < _numHierarchyLevels; ++i) {
            _hierarchies[i] = LinkedCells<Particle>(_boxMin, _boxMax, _cutoffs[i], skinPerTimestep, rebuildFrequency, cellSizeFactor, loadEstimator);
          }
        }
  
  // TODO: fix, add new option?
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::linkedCells; }

  [[nodiscard]] CellType getParticleCellTypeEnum() const override { return CellType::FullParticleCell; }

  // TODO: add reserve somehow
  void reserve(size_t numParticles, size_t numParticlesHaloEstimate) override {
    // do not do anything, because numParticles include all particles over all hierarchies
    // can do numParticles / _numHierarchyLevels, will not be accurate
  }

  void addParticleImpl(const ParticleType &p) override {
    _hierarchies[getHierarchyLevel(p)].addParticle(p);
  }

  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    _hierarchies[getHierarchyLevel(haloParticle)].addParticle(haloParticle);
  }

  // TODO: implement
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    // auto cells = _cellBlock.getNearbyHaloCells(haloParticle.getR(), this->getVerletSkin());
    // for (auto cellptr : cells) {
    //   bool updated = internal::checkParticleInCellAndUpdateByID(*cellptr, haloParticle);
    //   if (updated) {
    //     return true;
    //   }
    // }
    AutoPasLog(TRACE, "UpdateHaloParticle was not able to update particle: {}", haloParticle.toString());
    return false;
  }

  void deleteHaloParticles() override {
    for (const auto &linkedCells: _hierarchies) {
      linkedCells.deleteHaloParticles();
    }
  }

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _baseCutoff;
  double _skin;
  size_t _numHierarchyLevels;

  /**
   *
   * @param p Particle to add into HierarchialGrid
   * @return which Hierarchy the particle belongs to
   */
  size_t getHierarchyLevel(const ParticleType &p) {
    // scale size by baseCutoff
    const double cutoff = p.getSize() * _baseCutoff;
    for (size_t i = 0; i < _numHierarchyLevels; ++i) {
      if (_cutoffs[i] >= cutoff) {
        return i;
      }
    }
    AutoPasLog(ERROR, "Size of Particle times baseCutoff is bigger than biggest cutoff of HierarchicalGrid, "
                      "will result in wrong interaction calculation");
    return _numHierarchyLevels - 1;
  }

 protected:
  /**
   * load estimation algorithm for balanced traversals.
   */
  autopas::LoadEstimatorOption _loadEstimator;
  std::vector<autopas::LinkedCells<Particle>> _hierarchies;
  std::vector<double> _cutoffs;

};

}