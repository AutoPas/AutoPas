/**
 * @file HGC08SingleLevelTraversal.h
 * @author Alexander Glas
 * @date 12.04.2026
 */

#pragma once

#include "HGC08CellHandler.h"
#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCTraversalInterface.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the hg_c08 traversal, which is an adaptation of lc_c08 for Hierarchical Grids.
 * It calculates the intra-level interactions of one level, and the inter-level interactions of this level, with all
 * lower levels respectively. (e.g. when applied to level 3, it calculates 3-3, 3-2 and 3-1 interactions).
 *
 * The traversal uses the c08 base step performed on every single cell.
 * \image html C08.png "C08 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eight colors is applied.
 * \image html C08_domain.png "C08 domain coloring in 2D. 4 colors are required."
 * For every cell pair of the standard c08 base step, additionally to the intra-level interactions, each cell gets
 * decomposed into the lower level cells inside the domain of the higher level cell, and the inter-level interactions
 * are calculated.
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor>
class HGC08SingleLevelTraversal : public C08BasedTraversal<ParticleCell_T, PairwiseFunctor>,
                                  public LCTraversalInterface {
 public:
  /**
   * Constructor of the hg_c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin) of the higher level.
   * @param cellLength cell length of the higher level.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @param cellBlocks the CellBlocks of all levels, needed for inter-level interactions
   * @param interactionLengthsSquared the squared interaction lengths between the upperl level and all lower levels
   * respectively.
   * @param upperLevel The higher level.
   * @param haloRegionWidth the width of the halo region in length units, used to determine the number of halo cells to
   * be skipped in the traversal.
   * @param fittedGrids Whether the grids of the hierarchical grid are fitted to each other.
   */
  explicit HGC08SingleLevelTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                     double interactionLength, const std::array<double, 3> &cellLength,
                                     DataLayoutOption dataLayout, bool useNewton3,
                                     const std::vector<internal::CellBlock3D<ParticleCell_T> *> &cellBlocks,
                                     const std::vector<double> &interactionLengthsSquared, const size_t upperLevel,
                                     double haloRegionWidth, bool fittedGrids = true)
      : C08BasedTraversal<ParticleCell_T, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                           dataLayout, useNewton3),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap,
                     dataLayout, useNewton3, cellBlocks, interactionLengthsSquared, upperLevel, fittedGrids) {
    for (size_t d = 0; d < 3; ++d) {
      _haloRegionCellWidth[d] = std::ceil(haloRegionWidth / cellLength[d]);
    }
  }

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgridFit_c08; }

  /**
   * This traversal only makes sense for hierarchical grids, and only works if each cell size is at least twice as big
   * as the next lower level's cell size. However, this class is only used from within other traversals and only in
   * cases where it makes sense, so we always return true.
   * @return true
   */
  [[nodiscard]] bool isApplicable() const override { return true; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

 private:
  /*
   * The cell handler for managing the cells in the hierarchical grid.
   */
  HGC08CellHandler<ParticleCell_T, PairwiseFunctor> _cellHandler;
  /**
   * The width of the halo region in cells.
   */
  std::array<size_t, 3> _haloRegionCellWidth;
};

template <class ParticleCell_T, class PairwiseFunctor>
inline void HGC08SingleLevelTraversal<ParticleCell_T, PairwiseFunctor>::traverseParticles() {
  using namespace autopas::utils::ArrayMath::literals;

  const auto end = this->_cellsPerDimension - this->_haloRegionCellWidth;
  const auto stride = this->_overlap + 1ul;
  const auto offset = this->_haloRegionCellWidth - this->_overlap;

  this->colorTraversal(
      [&](unsigned long x, unsigned long y, unsigned long z) {
        const auto baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
        _cellHandler.processBaseCell(baseIndex);
      },
      end, stride, offset);
}
}  // namespace autopas
