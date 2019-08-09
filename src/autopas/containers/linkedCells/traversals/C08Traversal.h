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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class C08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
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
                        const double interactionLength, const std::array<double, 3> &cellLength)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

  void traverseParticlePairs() override;

  TraversalOption getTraversalType() const override { return TraversalOption::c08; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

  /**
   * C08 traversals are always usable.
   * @return
   */
  bool isApplicable() const override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    if (dataLayout == DataLayoutOption::cuda)
      return nDevices > 0;
    else
      return true;
  }

 private:
  C08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C08Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseParticlePairs() {
  auto &cells = *(this->_cells);
  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _cellHandler.processBaseCell(cells, baseIndex);
  });
}

}  // namespace autopas
