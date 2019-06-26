/**
 * @file AdaptiveLinkedCells.h
 *
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <autopas/cells/SortedCellView.h>
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/adaptiveLinkedCells/OctreeExternalNode.h"
#include "autopas/containers/adaptiveLinkedCells/OctreeInternalNode.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/linkedCells/traversals/LinkedCellTraversalInterface.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"
#include "autopas/options/DataLayoutOption.h"
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
class AdaptiveLinkedCells : public ParticleContainer<Particle, ParticleCell, SoAArraysType> {
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
                      const double cellSizeFactor = 1.0)
      : ParticleContainer<Particle, ParticleCell, SoAArraysType>(boxMin, boxMax, cutoff) {
    // set global values for tree
    OctreeExternalNode<Particle, ParticleCell>::setMaxElements(64);
    OctreeInternalNode<Particle, ParticleCell>::setMinElements(32);
    // compute cell length
    unsigned long _numCells = 1ul;
    for (int d = 0; d < 3; ++d) {
      const double diff = this->getBoxMax()[d] - this->getBoxMin()[d];
      auto cellsPerDim = static_cast<size_t>(std::floor(diff / (this->getCutoff() * cellSizeFactor)));
      // at least one central cell
      cellsPerDim = std::max(cellsPerDim, 1ul);

      //_cellLength[d] = diff / cellsPerDim;

      //_cellLengthReciprocal[d] = cellsPerDim / diff;  // compute with least rounding possible

      _numCells *= cellsPerDim;

      AutoPasLog(debug, "CellsPerDimension[{}]={}", d, cellsPerDim);
    }
    // +1, since last cell represents halo
    this->_cells.resize(_numCells + 1ul);
    // setup tree
    octree = std::make_unique<OctreeExternalNode<Particle, ParticleCell>>(
        this->_cells, 0, ArrayMath::mulScalar(ArrayMath::sub(boxMax, boxMin), 0.5), 0);
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

  void deleteHaloParticles() override { this->_cells.back().clear(); }

  ContainerOption getContainerType() override { return ContainerOption::adaptiveLinkedCells; }

  void addParticle(Particle &p) override {
    bool inBox = autopas::utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax());
    if (inBox) {
      ParticleCell &cell = octree->getContainingCell(p.getR());
      cell.addParticle(p);
    } else {
      utils::ExceptionHandler::exception(
          "LinkedCells: Trying to add a particle that is not inside the bounding box.\n" + p.toString());
    }
  }

  /**
   * @copydoc DirectSum::iteratePairwise
   */
  template <class ParticleFunctor, class Traversal>
  void iteratePairwise(ParticleFunctor *f, Traversal *traversal, bool useNewton3 = false) {
    AutoPasLog(debug, "Using traversal {}.", utils::StringUtils::to_string(traversal->getTraversalType()));

    traversal->initTraversal(this->_cells);
    if (auto *traversalInterface = dynamic_cast<LinkedCellTraversalInterface<ParticleCell> *>(traversal)) {
      traversalInterface->traverseCellPairs(this->_cells);

    } else {
      autopas::utils::ExceptionHandler::exception(
          "Trying to use a traversal of wrong type in AdaptiveLinkedCells::iteratePairwise. TraversalID: {}",
          traversal->getTraversalType());
    }
    traversal->endTraversal(this->_cells);
  }

  void updateContainer() override { /*octree = octree.update(); */
  }

  bool isContainerUpdateNeeded() override { return octree->isUpdateNeeded(); }

  TraversalSelectorInfo<ParticleCell> getTraversalSelectorInfo() override {
    return TraversalSelectorInfo<ParticleCell>({0ul, 0ul, 0ul}, this->getCutoff(), {0.0, 0.0, 0.0});
  }

  ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new internal::ParticleIterator<Particle, ParticleCell>(&this->_cells, 0, nullptr, behavior));
  }

  ParticleIteratorWrapper<Particle> getRegionIterator(const std::array<double, 3> &lowerCorner,
                                                      const std::array<double, 3> &higherCorner,
                                                      IteratorBehavior behavior = IteratorBehavior::haloAndOwned,
                                                      bool incSearchRegion = false) override {
    // @todo add implementation
    return ParticleIteratorWrapper<Particle>();
  }

  /**
   * returns reference to the data of AdaptiveLinkedCells
   * @return the data
   */
  std::vector<ParticleCell> &getCells() {
    // @todo add implementation
    return std::vector<ParticleCell>{};
  }

 protected:
  std::unique_ptr<OctreeNode<Particle, ParticleCell>> octree;
};

}  // namespace autopas
