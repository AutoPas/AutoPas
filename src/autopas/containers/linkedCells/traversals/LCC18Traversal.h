/**
 * @file LCC18Traversal.h
 * @author nguyen
 * @date 06.09.2018
 */

#pragma once

#include "LCTraversalInterface.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/containers/cellTraversals/C18BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the lc_c18 traversal.
 *
 * The traversal uses the c18 base step performed on every single cell.
 * \image html C18.png "C18 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eighteen colors is applied.
 * \image html C18_domain.png "C18 domain coloring in 2D. 6 colors are required."
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class LCC18Traversal : public C18BasedTraversal<ParticleCell, PairwiseFunctor>, public LCTraversalInterface {
 public:
  /**
   * Constructor of the lc_c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   * @todo Pass cutoff to _cellFunctor instead of interactionLength, unless this functor is used to build verlet-lists,
   * in that case the interactionLength is needed!
   */
  explicit LCC18Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                          const double interactionLength, const std::array<double, 3> &cellLength,
                          DataLayoutOption dataLayout, bool useNewton3)
      : C18BasedTraversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                         dataLayout, useNewton3),
        _cellFunctor(pairwiseFunctor, interactionLength /*should use cutoff here, if not used to build verlet-lists*/,
                     dataLayout, useNewton3) {
    this->computeOffsets();
  }

  void traverseParticles() override;

  /**
   * Computes all interactions between the base
   * cell and adjacent cells with greater a ID.
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c18; }

  /**
   * C18 traversal is always usable.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override { return true; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellFunctor.setSortingThreshold(sortingThreshold); }

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  internal::CellFunctor<ParticleCell, PairwiseFunctor,
                        /*bidirectional*/ true>
      _cellFunctor;
};

template <class ParticleCell, class PairwiseFunctor>
void LCC18Traversal<ParticleCell, PairwiseFunctor>::processBaseCell(std::vector<ParticleCell> &cells, unsigned long x,
                                                                    unsigned long y, unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  const unsigned long xArray = this->getIndex(x, 0);
  const unsigned long yArray = this->getIndex(y, 1);

  ParticleCell &baseCell = cells[baseIndex];
  auto &offsets = this->_cellOffsets[yArray][xArray];
  for (auto const &[offset, r] : offsets) {
    unsigned long otherIndex = baseIndex + offset;
    ParticleCell &otherCell = cells[otherIndex];

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseCell);
    } else {
      this->_cellFunctor.processCellPair(baseCell, otherCell, r);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor>
inline void LCC18Traversal<ParticleCell, PairwiseFunctor>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->template c18Traversal</*allCells*/ false>(
      [&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
