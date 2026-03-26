/**
 * @file LCC04CombinedSoATraversal.h
 * @author C. Menges
 * @date 02.06.2019
 */

#pragma once

#include "LCC04SoACellHandler.h"
#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/C04BasedTraversal.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c04 traversal.
 *
 * The traversal uses the c04 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eight colors is applied.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 */
template <class ParticleCell, class PairwiseFunctor>
class LCC04CombinedSoATraversal : public C04BasedTraversal<ParticleCell, PairwiseFunctor, 2>,
                                  public LCTraversalInterface {
 public:
  /**
   * Constructor of the c04 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length.
   * @param cellLength cell length.
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit LCC04CombinedSoATraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                     const double interactionLength, const std::array<double, 3> &cellLength,
                                     DataLayoutOption dataLayout, bool useNewton3)
      : C04BasedTraversal<ParticleCell, PairwiseFunctor, 2>(dims, pairwiseFunctor, interactionLength, cellLength,
                                                            dataLayout, useNewton3),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, dataLayout, useNewton3,
                     this->_overlap) {}

  void traverseParticles() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c04_combined_SoA; }

  /**
   * lc_c04_combined_SoA traversals are only usable with dataLayout SoA.
   * @todo Currently there is a bug when there is an overlap of more than one cell (typically due to CSF<1.0):
   * https://github.com/AutoPas/AutoPas/issues/354
   * once this bug is fixed, reenable this traversal again for arbitrary `_overlap`s.
   * @return
   */
  [[nodiscard]] bool isApplicableToDomain() const override {
    return this->_overlap[0] == 1 and this->_overlap[1] == 1 and this->_overlap[2] == 1;
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  LCC04SoACellHandler<ParticleCell, PairwiseFunctor> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor>
inline void LCC04CombinedSoATraversal<ParticleCell, PairwiseFunctor>::traverseParticles() {
  _cellHandler.resizeBuffers();
  auto &cells = *(this->_cells);
  this->c04Traversal(
      [&](unsigned long x, unsigned long y, unsigned long z) { _cellHandler.processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
