/**
 * @file C04Traversal.h
 * @author C. Menges
 * @date 02.06.2019
 */

#pragma once

#include "C04CellHandler.h"
#include "LinkedCellTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C04BasedTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
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
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C04Traversal : public C04BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, 2>,
                     public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c04 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit C04Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                        const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : C04BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3, 2>(dims, pairwiseFunctor, cutoff,
                                                                                    cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, cutoff, cellLength, this->_overlap) {}

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;
  TraversalOption getTraversalType() override { return TraversalOption::c04SoA; }

  /**
   * c04SoA traversals are only usable with DataLayout SoA.
   * @return
   */
  bool isApplicable() override { return DataLayout == DataLayoutOption::soa; }

 private:
  C04CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C04Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  _cellHandler.resizeBuffers();
  this->c04Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, x, y, z);
  });
}

}  // namespace autopas
