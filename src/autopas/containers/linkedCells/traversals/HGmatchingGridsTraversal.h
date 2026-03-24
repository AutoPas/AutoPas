/**
 * @file HGmatchingGridsTraversal.h
 * @author AlexanderGlas
 * @date .02.2026
 */

#pragma once

#include "HGC08Traversal.h"
#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {
/**
 * @todo own documentation
 * For each level, LCC08Traversal is used. For the cross-level interactions, for each level x only smaller levels
 * are iterated (newton3 on only). The cells on level x are iterated with colors (dynamic color count based on ratio
 * of cell lengths between level x and y) so that the cells on the lower level y
 * that are considered for each cell on level x do not intersect.
 * To reduce number of colors and increase memory efficiency, instead of only 1 upper
 * level cell a block of cells is assigned to a thread at a time. The size of block is calculated dynamically
 * by considering upper and lower cell lengths and number of threads. The number of blocks per color is at least
 * num_threads * 4 or 8, depending on the option.
 * @tparam ParticleCell_T type of Particle cell
 * @tparam Functor_T type of Functor
 */
template <class ParticleCell_T, class Functor_T>
class HGmatchingGridsTraversal : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGmatchingGridsTraversal(Functor_T *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell_T>(dataLayout, useNewton3),
        _functor(functor),
        _dataLayoutConverter(_functor, dataLayout) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported");
    }

    const double haloRegionLength = this->_cutoffs.back() + this->_skin;

    // get a vector of cell blocks for easier access in the traversal
    std::vector<CellBlock *> cellBlocks;
    cellBlocks.reserve(this->_levels->size());
    for (auto &linkedCell : *this->_levels) {
      cellBlocks.emplace_back(&linkedCell->getCellBlock());
    }

    // computeInteractions across different levels
    for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
      // InteractionLengths squared from upper toall lower levels, used for checking to skip cells out of range
      // @todo see at what point that is worth
      std::vector<double> interactionLengthsSquared(upperLevel);
      for (size_t lowerLevel = 0; lowerLevel < upperLevel; lowerLevel++) {
        auto tmp = this->getInteractionLength(upperLevel, lowerLevel);
        interactionLengthsSquared[lowerLevel] = tmp * tmp;
      }

      // prepare and perform HGc08 traversal
      // set actual halo region length of the traversal
      // this will let traversal skip going through unnecessary halo cells according to atacan
      auto currentTraversal = this->generateNewMixedTraversal(upperLevel, cellBlocks, interactionLengthsSquared);
      currentTraversal->setHaloRegionLength(haloRegionLength);
      currentTraversal->traverseParticles();
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_matching; };

  // only if matching
  [[nodiscard]] bool isApplicable() const override { return true; }

  /**
   * load Data Layouts required for this Traversal if cells have been set through setCellsToTraverse().
   */
  void initTraversal() override {
    for (auto &level : *this->_levels) {
      auto &cells = level->getCells();
      AUTOPAS_OPENMP(parallel for schedule(dynamic, 1))
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.loadDataLayout(cells[i]);
      }
    }
  }

  /**
   * write Data to AoS if cells have been set through setCellsToTraverse().
   */
  void endTraversal() override {
    for (auto &level : *this->_levels) {
      auto &cells = level->getCells();
      AUTOPAS_OPENMP(parallel for schedule(dynamic, 1))
      for (size_t i = 0; i < cells.size(); ++i) {
        _dataLayoutConverter.storeDataLayout(cells[i]);
      }
    }
  }

 protected:
  // HGTraversalBase for some reason needs hardcoded FullParticleCell
  using CellBlock = internal::CellBlock3D<FullParticleCell<Particle>>;

  Functor_T *_functor;
  utils::DataLayoutConverter<Functor_T> _dataLayoutConverter;
  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  // TODO: haloLength vs interaction length
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) override {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    return std::make_unique<LCC08Traversal<ParticleCell_T, Functor_T>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }

  std::unique_ptr<HGC08Traversal<FullParticleCell<Particle>, Functor_T>> generateNewMixedTraversal(
      const size_t upper, const std::vector<CellBlock *> &cellBlocks,
      const std::vector<double> &interactionLengthsSquared) {
    const auto traversalInfo = this->getTraversalSelectorInfo(upper);
    return std::make_unique<HGC08Traversal<FullParticleCell<Particle>, Functor_T>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3, cellBlocks, interactionLengthsSquared, upper, this->_fittedGrids);
  }
};
}  // namespace autopas