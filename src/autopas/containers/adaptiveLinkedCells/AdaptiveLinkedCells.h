/**
 * @file AdaptiveLinkedCells.h
 *
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <autopas/cells/SortedCellView.h>
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/adaptiveLinkedCells/Octree.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LinkedCellTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/selectors/TraversalSelectorInfoAdaptive.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StringUtils.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * AdaptiveLinkedCells class.
 * This class uses a list of neighboring cells to store the particles.
 * These cells dimensions are at least as large as the given cutoff radius,
 * therefore short-range interactions only need to be calculated between
 * particles in neighboring cells.
 * All cells are managed by an octree.
 * @tparam Particle type of the particles that need to be stored
 * @tparam ParticleCell type of the ParticleCells that are used to store the particles
 * @tparam SoAArraysType type of the SoA, needed for verlet lists
 */
template <class Particle, class ParticleCell, class SoAArraysType = typename Particle::SoAArraysType>
class AdaptiveLinkedCells final : public ParticleContainer<Particle, ParticleCell, SoAArraysType> {
 public:
  /**
   * Constructor of the AdaptiveLinkedCells class.
   * @param boxMin
   * @param boxMax
   * @param cutoff
   * @param cellSizeFactor Cell size factor relative to cutoff.
   * By default all applicable traversals are allowed.
   */
  AdaptiveLinkedCells(const std::array<double, 3> boxMin, const std::array<double, 3> boxMax, const double cutoff,
                      const double skin, const double cellSizeFactor = 1.0)
      : ParticleContainer<Particle, ParticleCell, SoAArraysType>(boxMin, boxMax, cutoff, skin),
        octree(this->_cells, boxMin, boxMax) {
    // set global values for tree
    Octree<Particle, ParticleCell>::setMaxElements(64);
    Octree<Particle, ParticleCell>::setMinElements(32);
    // compute cell length

    for (unsigned int d = 0; d < 3; ++d) {
      _cellLength[d] = this->getBoxMax()[d] - this->getBoxMin()[d];
      _cellsPerDimension[d] = static_cast<size_t>(std::floor(_cellLength[d] / (this->getCutoff() * cellSizeFactor)));
      // at least one central cell
      _cellsPerDimension[d] = std::max(_cellsPerDimension[d], 1ul);
      // round down to next exponential of 2
      _cellsPerDimension[d] = 1 << static_cast<int>(log2(_cellsPerDimension[d]));
    }
    // make sure all axis have the same size
    _cellsPerDimension.fill(*std::min_element(_cellsPerDimension.cbegin(), _cellsPerDimension.cend()));
    AutoPasLog(debug, "CellsPerDimension[{}, {}, {}]", _cellsPerDimension[0], _cellsPerDimension[1],
               _cellsPerDimension[2]);

    for (unsigned int d = 0; d < 3; ++d) {
      _cellLength[d] /= _cellsPerDimension[d];
    }

    size_t numCells = _cellsPerDimension[0] * _cellsPerDimension[1] * _cellsPerDimension[2];
    // +1, since last cell represents halo
    this->_cells.resize(numCells + 1ul);
    octree.init(_cellsPerDimension);
  }

  void addHaloParticle(Particle &haloParticle) override {
    bool inHalo = not utils::inBox(haloParticle.getR(), this->getBoxMin(), this->getBoxMax());
    if (inHalo) {
      ParticleCell &cell = this->_cells.back();
      cell.addParticle(haloParticle);
    } else {
      utils::ExceptionHandler::exception(
          "LinkedCells: Trying to add a halo particle that is actually in the bounding box.\n" +
          haloParticle.toString());
    }
  }

  bool updateHaloParticle(Particle &haloParticle) override {
    /// @todo add impl
    return false;
  }
  void rebuildNeighborLists(TraversalInterface *traversal) override {
    // nothing to do.
  }

  void deleteHaloParticles() override { this->_cells.back().clear(); }

  ContainerOption getContainerType() override { return ContainerOption::adaptiveLinkedCells; }

  void addParticle(Particle &p) override {
    bool inBox = autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      ParticleCell &cell = octree.getContainingCell(p.getR());
      cell.addParticle(p);
    } else {
      utils::ExceptionHandler::exception(
          "LinkedCells: Trying to add a particle that is not inside the bounding box.\n" + p.toString());
    }
    AutoPasLog(debug, "Number of particles {} ", getCells()[0].numParticles());
  }

  /**
   * @copydoc DirectSum::iteratePairwise
   */
  void iteratePairwise(TraversalInterface *traversal) override {
    AutoPasLog(debug, "Using traversal {}.", utils::StringUtils::to_string(traversal->getTraversalType()));

    traversal->initTraversal();
    /*if (auto *traversalInterface = dynamic_cast<LinkedCellTraversalInterface<ParticleCell> *>(traversal)) {
      traversalInterface->traverseCellPairs(this->_cells);

    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in AdaptiveLinkedCells::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }*/
    traversal->endTraversal();
  }

  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<Particle> updateContainer() override {
    octree.update();
    return std::vector<Particle>{};
  }

  bool isContainerUpdateNeeded() override { return octree.isUpdateNeeded(); }

  std::unique_ptr<TraversalSelectorInfo> getTraversalSelectorInfo() override {
    return std::make_unique<TraversalSelectorInfoAdaptive>(_cellsPerDimension, this->getCutoff(),
                                                           _cellLength /*, octree*/);
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_cells, 0, nullptr, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    // @todo add implementation
    return ParticleIteratorWrapper<Particle>();
  }

  /**
   * returns reference to the data of AdaptiveLinkedCells
   * @return the data
   */
  std::vector<ParticleCell> &getCells() { return this->_cells; }

 protected:
  Octree<Particle, ParticleCell> octree;
  std::array<unsigned long, 3> _cellsPerDimension;
  std::array<double, 3> _cellLength;
};

}  // namespace autopas
