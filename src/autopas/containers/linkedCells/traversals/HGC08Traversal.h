/**
 * @file HGC08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "HGC08CellHandler.h"
#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the lc_c08 traversal.
 *
 * The traversal uses the c08 base step performed on every single cell.
 * \image html C08.png "C08 base step in 2D. (dark blue cell = base cell)"
 * Since these steps overlap a domain coloring with eight colors is applied.
 * \image html C08_domain.png "C08 domain coloring in 2D. 4 colors are required."
 *
 * @tparam ParticleCell_T the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell_T, class PairwiseFunctor>
class HGC08Traversal : public C08BasedTraversal<ParticleCell_T, PairwiseFunctor>, public LCTraversalInterface {
 public:
  /**
   * Constructor of the lc_c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit HGC08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                          double interactionLength, const std::array<double, 3> &cellLength,
                          DataLayoutOption dataLayout, bool useNewton3,
                          const std::vector<internal::CellBlock3D<ParticleCell_T> *> &cellBlocks,
             const std::vector<double> &interactionLengthsSquared, const size_t upperLevel,
             bool fittedGrids = false)
      : C08BasedTraversal<ParticleCell_T, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                           dataLayout, useNewton3),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap,
           dataLayout, useNewton3, cellBlocks, interactionLengthsSquared, upperLevel, fittedGrids) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_matching; }

  /**
   * C08 traversals are always usable.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override { return true; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

 private:
  HGC08CellHandler<ParticleCell_T, PairwiseFunctor> _cellHandler;
};

template <class ParticleCell_T, class PairwiseFunctor>
inline void HGC08Traversal<ParticleCell_T, PairwiseFunctor>::traverseParticles() {
  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(baseIndex);
  });
}

}  // namespace autopas
