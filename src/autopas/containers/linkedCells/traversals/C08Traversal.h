/**
 * @file C08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "C08CellHandler.h"
#include "LinkedCellTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c08 traversal.
 *
 * The traversal uses the c08 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eight colors is applied.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>,
                     public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   */
  explicit C08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                        const double interactionLength = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;
  TraversalOption getTraversalType() override { return TraversalOption::c08; }

  /**
   * C08 traversals are always usable.
   * @return
   */
  bool isApplicable() override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    if (DataLayout == DataLayoutOption::cuda)
      return nDevices > 0;
    else
      return true;
  }

 private:
  C08CellHandler<ParticleCell, PairwiseFunctor, DataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, baseIndex);
  });
}

}  // namespace autopas
