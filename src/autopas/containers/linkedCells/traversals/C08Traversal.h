/**
 * @file C08Traversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "C08LikeBaseCellProcessor.h"
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
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C08Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>,
                     public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C08Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>(dims, pairwiseFunctor),
        _baseCellProcessor(pairwiseFunctor, this->_cellsPerDimension) {}

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;
  TraversalOptions getTraversalType() override { return TraversalOptions::c08; }

 private:
  C08LikeBaseCellProcessor<ParticleCell, PairwiseFunctor, useSoA, useNewton3> _baseCellProcessor;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C08Traversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
    unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
    _baseCellProcessor.processBaseCell(cells, baseIndex);
  });
}

}  // namespace autopas