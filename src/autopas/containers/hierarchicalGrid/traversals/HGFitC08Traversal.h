/**
 * @file HGFitC08Traversal.h
 * @author Alexander Glas
 * @date 12.04.2026
 */

#pragma once

#include "HGC08SingleLevelTraversal.h"
#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {
/**
 * This class provides the hgridFit_c08 traversal.
 *
 * The traversal calculates the pairwise interactions of a Hierarchical Grid, by iterating over all levels and applying
 * hgc08SingleLevelTraversal, which calculates the intra-level interactions for one specific level and the inter-level
 * interactions between this level and all lower levels respectively.
 *
 * @tparam ParticleCell_T type of Particle cell
 * @tparam Functor_T type of Functor
 */
template <class ParticleCell_T, class Functor_T>
class HGFitC08Traversal : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  /**
   * Constructor of the lc_c08 traversal.
   * @param functor The functor that defines the interaction of two particles.
   * @param numLevels Number of levels in the hierarchical grid.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit HGFitC08Traversal(Functor_T *functor, size_t numLevels, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell_T>(numLevels, dataLayout, useNewton3),
        _functor(functor),
        _dataLayoutConverter(_functor, dataLayout) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgridFit_c08; };

  /**
   * C08 traversals are always usable for fitted grids.
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
  using CellBlock = internal::CellBlock3D<FullParticleCell<Particle>>;
  Functor_T *_functor;
  utils::DataLayoutConverter<Functor_T> _dataLayoutConverter;

  /**
   * HGFitc08Traversal does intra and inter-level traversals combined, so generating seperate intra-level traversals is
   * not needed.
   * @param level which HGrid level to generate a traversal for
   * @return
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) override { return nullptr; }
};

template <class ParticleCell_T, class Functor_T>
inline void HGFitC08Traversal<ParticleCell_T, Functor_T>::traverseParticles() {
  const double haloRegionWidth = this->_maxCutoffPerLevel.back() + this->_skin;

  // get a vector of cell blocks for easier access in the traversal
  std::vector<CellBlock *> cellBlocks;
  cellBlocks.reserve(this->_levels->size());
  for (auto &linkedCell : *this->_levels) {
    cellBlocks.emplace_back(&linkedCell->getCellBlock());
  }

  // computeInteractions across different levels
  for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
    // InteractionLengths squared from upper to all lower levels, used to skip cells out of range
    std::vector<double> interactionLengthsSquared(upperLevel);
    for (size_t lowerLevel = 0; lowerLevel < upperLevel; lowerLevel++) {
      auto interactionLength = this->getInteractionLength(upperLevel, lowerLevel);
      interactionLengthsSquared[lowerLevel] = interactionLength * interactionLength;
    }

    // prepare and perform HGc08 traversal
    const auto traversalInfo = this->getTraversalSelectorInfo(upperLevel);
    auto currentTraversal = std::make_unique<HGC08SingleLevelTraversal<FullParticleCell<Particle>, Functor_T>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3, cellBlocks, interactionLengthsSquared, upperLevel, haloRegionWidth);
    currentTraversal->traverseParticles();
  }
}

}  // namespace autopas