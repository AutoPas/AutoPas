/**
 * @file LCC08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "LCC08CellHandler.h"
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
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class LCC08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor>, public LCTraversalInterface {
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
  explicit LCC08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                          double interactionLength, const std::array<double, 3> &cellLength,
                          DataLayoutOption dataLayout, bool useNewton3)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                         dataLayout, useNewton3),

        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap,
                     dataLayout, useNewton3) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c08; }

  /**
   * LC C08 is always applicable to the domain.
   * @return true
   */
  [[nodiscard]] bool isApplicableToDomain() const override { return true; }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   */
  void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

 private:
  LCC08CellHandler<ParticleCell, PairwiseFunctor> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor>
inline void LCC08Traversal<ParticleCell, PairwiseFunctor>::traverseParticles() {
  auto &cells = *(this->_cells);
  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, baseIndex);
  });
}

}  // namespace autopas
