/**
 * @file VerletListsCells.h
 * @author nguyen
 * @date 30.08.18
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include "autopas/containers/LoadEstimators.h"
#include "autopas/containers/cellPairTraversals/BalancedTraversal.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletListsCellBased/VerletListsLinkedBase.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCellsHelpers.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/selectors/TraversalSelector.h"

namespace autopas {

/**
 * Linked Cells with Verlet Lists container.
 * The VerletListsCells class uses neighborhood lists for each cell
 * to calculate pairwise interactions of particles.
 * It is optimized for a constant, i.e. particle independent, cutoff radius of
 * the interaction.
 * Cells are created using a cell size of at least cutoff + skinPerTimestep*rebuildFrequency.
 * @tparam Particle
 * @tparam NeighborList The neighbor list used by this container.
 */

template <class Particle, class NeighborList>
class VerletListsCells : public VerletListsLinkedBase<Particle> {
  using verlet_internal = VerletListsCellsHelpers<FullParticleCell<Particle>>;
  using ParticleCell = FullParticleCell<Particle>;

 public:
  /**
   * Constructor of the VerletListsCells class.
   * The neighbor lists are build using a search radius of cutoff + skin*rebuildfrequency.
   * @param boxMin the lower corner of the domain
   * @param boxMax the upper corner of the domain
   * @param cutoff the cutoff radius of the interaction
   * @param rebuildFrequency the rebuild Frequency
   * @param skinPerTimestep the skin radius per Timestep
   * @param cellSizeFactor cell size factor relative to cutoff
   * @param loadEstimator load estimation algorithm for balanced traversals
   * @param buildType data layout of the particles which are used to generate the neighbor lists
   */
  VerletListsCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                   const double skinPerTimestep = 0, const unsigned int rebuildFrequency = 2,
                   const double cellSizeFactor = 1.0,
                   const LoadEstimatorOption loadEstimator = LoadEstimatorOption::squaredParticlesPerCell,
                   typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value buildType =
                       VerletListsCellsHelpers<Particle>::VLCBuildType::soaBuild)
      : VerletListsLinkedBase<Particle>(boxMin, boxMax, cutoff, skinPerTimestep, rebuildFrequency,
                                        compatibleTraversals::allVLCCompatibleTraversals(), cellSizeFactor),
        _loadEstimator(loadEstimator),
        _buildType(buildType) {}

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] virtual ContainerOption getContainerType() const { return _neighborList.getContainerType(); }

  /**
   * Generates the load estimation function depending on _loadEstimator.
   * @return load estimator function object.
   */
  BalancedTraversal::EstimatorFunction getLoadEstimatorFunction() {
    switch (this->_loadEstimator) {
      case LoadEstimatorOption::squaredParticlesPerCell: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          return loadEstimators::squaredParticlesPerCell((this->_linkedCells).getCells(), cellsPerDimension,
                                                         lowerCorner, upperCorner);
        };
      }
      case LoadEstimatorOption::neighborListLength: {
        return [&](const std::array<unsigned long, 3> &cellsPerDimension,
                   const std::array<unsigned long, 3> &lowerCorner, const std::array<unsigned long, 3> &upperCorner) {
          return loadEstimators::neighborListLength<Particle, NeighborList>(_neighborList, cellsPerDimension,
                                                                            lowerCorner, upperCorner);
        };
      }

      case LoadEstimatorOption::none:
        [[fallthrough]];
      default: {
        return
            [&](const std::array<unsigned long, 3> &cellsPerDimension, const std::array<unsigned long, 3> &lowerCorner,
                const std::array<unsigned long, 3> &upperCorner) { return 1; };
      }
    }
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    // Check if traversal is allowed for this container and give it the data it needs.
    _neighborList.setUpTraversal(traversal);
    if (auto *balancedTraversal = dynamic_cast<BalancedTraversal *>(traversal)) {
      balancedTraversal->setLoadEstimator(getLoadEstimatorFunction());
    }

    /*
    // count number of current dummy particles
    size_t dummyCounter {0};
    for (auto pIter = this->begin(); pIter != this->end(); ++pIter) {
      if (pIter->isDummy()) {
        ++dummyCounter;
      }
    }
    _dummyParticleCounter.push_back(dummyCounter);
    */

    traversal->initTraversal();
    traversal->traverseParticlePairs();
    traversal->endTraversal();
  }

  /**
   * Gets the number of neighbors over all neighbor lists that belong to this particle.
   * @param particle
   * @return the size of the neighbor list(s) of this particle
   */
  size_t getNumberOfPartners(const Particle *particle) const { return _neighborList.getNumberOfPartners(particle); }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    this->_verletBuiltNewton3 = traversal->getUseNewton3();

    _neighborList.buildAoSNeighborList(this->_linkedCells, this->_verletBuiltNewton3, this->getCutoff(),
                                       this->getVerletSkin(), this->getInteractionLength(), TraversalOption::lc_c01_b08,
                                       _buildType, this->_partialRebuilding);

    if (traversal->getDataLayout() == DataLayoutOption::soa) {
      _neighborList.generateSoAFromAoS(this->_linkedCells);
    }

    if (!this->_partialRebuilding) {
      this->_cellRebuildCounter.emplace_back(this->_linkedCells.getCells().size());
    }
    // the neighbor list is now valid
    this->_neighborListIsValid.store(true, std::memory_order_relaxed);
  }

  virtual bool neighborListsAreValid() override { return true; }

  /**
   * Return the cell length of the underlying linked cells structure, normally needed only for unit tests.
   * @return
   */
  [[nodiscard]] const std::array<double, 3> &getCellLength() const {
    return this->_linkedCells.getCellBlock().getCellLength();
  }

  size_t getNumRebuildCells() override {
    return std::accumulate(this->_cellRebuildCounter.begin(), this->_cellRebuildCounter.end(), 0);
  }

  size_t getNumInflowCells() override {
    return std::accumulate(this->_cellInflowCounter.begin(), this->_cellInflowCounter.end(), 0);
  }

  size_t getNumOutflowCells() override {
    return std::accumulate(this->_cellOutflowCounter.begin(), this->_cellOutflowCounter.end(), 0);
  }

  size_t getNumMovingCells() override {
    return std::accumulate(this->_cellMovingCounter.begin(), this->_cellMovingCounter.end(), 0);
  }

  size_t getNumDummyParticles() override {
    return std::accumulate(this->_dummyParticleCounter.begin(), this->_dummyParticleCounter.end(), 0);
  }

 protected:
  /**
   * Neighbor list abstraction for neighbor list used in the container.
   */
  NeighborList _neighborList;

  std::vector<size_t> _cellRebuildCounter;

  std::vector<size_t> _cellInflowCounter;

  std::vector<size_t> _cellOutflowCounter;

  std::vector<size_t> _cellMovingCounter;

  std::vector<size_t> _dummyParticleCounter;

private:
  /**
   * Load estimation algorithm for balanced traversals.
   */
  autopas::LoadEstimatorOption _loadEstimator;

  /**
   * Data layout of the particles which are used to generate the neighbor lists.
   */
  typename VerletListsCellsHelpers<Particle>::VLCBuildType::Value _buildType;
};
}  // namespace autopas
