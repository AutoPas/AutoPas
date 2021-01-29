/**
 * @file LCC08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "LCC08CellHandler.h"
#include "LCTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                       public LCTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the lc_c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit LCC08Traversal(const std::array<uint64_t, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                          const double interactionLength, const std::array<double, 3> &cellLength)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

  void traverseParticlePairs() override;

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c08; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  /**
   * C08 traversals are always usable.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override { return not(dataLayout == DataLayoutOption::cuda); }

 private:
  LCC08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  this->c08Traversal([&](uint64_t x, uint64_t y, uint64_t z) {
    uint64_t baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, baseIndex);
  });
}

}  // namespace autopas
