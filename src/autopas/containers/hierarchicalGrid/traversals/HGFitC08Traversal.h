/**
 * @file HGFitC08Traversal.h
 * @author Alexander Glas
 * @date 12.04.2026
 */

#pragma once

#include "HGC08CellHandler.h"
#include "HGTraversalBase.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {
/**
 * This class provides the hgridFit_c08 traversal.
 *
 * The traversal calculates the pairwise interactions of a Hierarchical Grid, by iterating over all levels calculating
 * the intra-level interactions for this specific level and the inter-level interactions between this level and all
 * lower levels respectively.
 *
 * This per-level traversal uses the c08 base step performed on every single cell.
 * \image html C08.png "C08 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eight colors is applied.
 * \image html C08_domain.png "C08 domain coloring in 2D. 4 colors are required."
 * For every cell pair of the standard c08 base step, additionally to the intra-level interactions, each cell gets
 * decomposed into the lower level cells inside the domain of the higher level cell, and the inter-level interactions
 * are calculated.
 * Also has an intra-level only mode, which works just like LCC08Traversal, but considers the potentially larger halo
 * regions of HierarchicalGrid levels to skip unnecessary halo-halo interactions.
 *
 * @tparam ParticleCell_T type of Particle cell
 * @tparam Functor_T type of Functor
 */
template <class ParticleCell_T, class Functor_T>
class HGFitC08Traversal : public HGTraversalBase<ParticleCell_T>,
                          public C08BasedTraversal<ParticleCell_T, Functor_T>,
                          public LCTraversalInterface {
 public:
  /**
   * Particle type handled by this traversal.
   */
  using Particle = typename ParticleCell_T::ParticleType;

  /**
   * Constructor of the hgridFit_c08 traversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param numLevels Number of levels in the hierarchical grid.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @param intraLevel Non-negative integer indicates to only calculate intra-level interactions for that specific
   * level. Works just like LCC08Traversal, but considers the potentially larger halo regions of HierarchicalGrid levels
   * to skip unnecessary halo-halo interactions.
   */
  explicit HGFitC08Traversal(Functor_T *pairwiseFunctor, size_t numLevels, DataLayoutOption dataLayout,
                             const bool useNewton3, const int intraLevel = -1)
      : HGTraversalBase<ParticleCell_T>(numLevels, dataLayout, useNewton3),
        C08BasedTraversal<ParticleCell_T, Functor_T>({}, pairwiseFunctor, 0, {}, dataLayout, useNewton3),
        _pairwiseFunctor(*pairwiseFunctor),
        _dataLayoutConverter(pairwiseFunctor, dataLayout),
        _intraLevel(intraLevel) {
    if (intraLevel >= static_cast<int>(numLevels)) {
      autopas::utils::ExceptionHandler::exception("HGFitC08Traversal: intraLevel must be smaller than numLevels.");
    }
  }

  void traverseParticles() override;

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _sortingThreshold = sortingThreshold; }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgridFit_c08; };

  /**
   * C08 traversals are always usable for fitted grids.
   * @todo: In the future sort out single level Hgrids.
   * @return true
   */
  [[nodiscard]] bool isApplicable() const override { return true; }

  /**
   * Load Data Layouts required for this Traversal.
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
   * Write Data to AoS.
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
  /**
   * Cell block type used for hierarchical-grid traversal.
   */
  using CellBlock = internal::CellBlock3D<FullParticleCell<Particle>>;
  /**
   * Pairwise functor used to compute interactions.
   */
  Functor_T &_pairwiseFunctor;
  /**
   * Utility that converts between AoS and SoA data layouts for the functor.
   */
  utils::DataLayoutConverter<Functor_T> _dataLayoutConverter;
  /**
   * Non-negative integer indicates to only perform intra-level interactions for that level.
   */
  const int _intraLevel;
  /**
   * Calculates the intra-level interactions for the level indicated by _intraLevel. Works just like LCC08Traversal, but
   * considers the potentially larger halo regions of HierarchicalGrid levels to skip unnecessary halo-halo
   * interactions.
   */
  void intraLevelOnlyTraversal();

  /**
   * Sorting threshold forwarded to local cell handlers.
   */
  size_t _sortingThreshold = 0;
};

template <class ParticleCell_T, class Functor_T>
inline void HGFitC08Traversal<ParticleCell_T, Functor_T>::traverseParticles() {
  using namespace autopas::utils::ArrayMath::literals;
  if (this->_intraLevel >= 0) {
    intraLevelOnlyTraversal();
    return;
  }
  // get a vector of cell blocks for easier access in the traversal
  std::vector<CellBlock *> cellBlocks;
  cellBlocks.reserve(this->_levels->size());
  for (auto &linkedCell : *this->_levels) {
    cellBlocks.emplace_back(&linkedCell->getCellBlock());
  }

  // compute top-down interactions for each level
  for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
    // InteractionLengths squared from upper to all lower levels, used to skip cells out of range
    std::vector<double> interactionLengthsSquared(upperLevel);
    for (size_t lowerLevel = 0; lowerLevel < upperLevel; lowerLevel++) {
      auto interactionLength = this->getInteractionLength(upperLevel, lowerLevel);
      interactionLengthsSquared[lowerLevel] = interactionLength * interactionLength;
    }
    // prepare for grid of current upperLevel
    this->changeGrid(this->_maxCutoffPerLevel[upperLevel] + this->_skin, cellBlocks[upperLevel]->getCellLength(),
                     cellBlocks[upperLevel]->getCellsPerDimensionWithHalo());
    HGC08CellHandler cellHandler{&_pairwiseFunctor,
                                 this->_cellsPerDimension,
                                 this->_interactionLength,
                                 this->_cellLength,
                                 this->_overlap,
                                 this->_dataLayout,
                                 this->_useNewton3,
                                 cellBlocks,
                                 interactionLengthsSquared,
                                 upperLevel};
    cellHandler.setSortingThreshold(_sortingThreshold);

    // Perform c08 based traversal, but consider the potentially larger halo regions of HierarchicalGrid levels to skip
    // unnecessary halo-halo interactions
    const unsigned long haloRegionCellWidth = cellBlocks[upperLevel]->getCellsPerInteractionLength();
    const auto end = this->_cellsPerDimension - haloRegionCellWidth;
    const auto stride = this->_overlap + 1ul;
    const std::array<unsigned long, 3> offset =
        std::array<unsigned long, 3>{haloRegionCellWidth, haloRegionCellWidth, haloRegionCellWidth} - this->_overlap;
    this->colorTraversal(
        [&](unsigned long x, unsigned long y, unsigned long z) {
          const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
          cellHandler.processBaseCell(baseIndex);
        },
        end, stride, offset);
  }
}

template <class ParticleCell_T, class Functor_T>
inline void HGFitC08Traversal<ParticleCell_T, Functor_T>::intraLevelOnlyTraversal() {
  using namespace autopas::utils::ArrayMath::literals;
  CellBlock &cellBlock = this->_levels->at(this->_intraLevel)->getCellBlock();
  this->changeGrid(this->_maxCutoffPerLevel[this->_intraLevel] + this->_skin, cellBlock.getCellLength(),
                   cellBlock.getCellsPerDimensionWithHalo());
  LCC08CellHandler<ParticleCell_T, Functor_T> cellHandler{
      &_pairwiseFunctor, this->_cellsPerDimension, this->_interactionLength, this->_cellLength,
      this->_overlap,    this->_dataLayout,        this->_useNewton3};
  cellHandler.setSortingThreshold(_sortingThreshold);
  const unsigned long haloRegionCellWidth = cellBlock.getCellsPerInteractionLength();
  const auto end = this->_cellsPerDimension - haloRegionCellWidth;
  const auto stride = this->_overlap + 1ul;
  const std::array<unsigned long, 3> offset =
      std::array<unsigned long, 3>{haloRegionCellWidth, haloRegionCellWidth, haloRegionCellWidth} - this->_overlap;
  this->colorTraversal(
      [&](unsigned long x, unsigned long y, unsigned long z) {
        const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
        cellHandler.processBaseCell(this->_levels->at(this->_intraLevel)->getCells(), baseIndex);
      },
      end, stride, offset);
}
}  // namespace autopas