/**
 * @file SlicedBlkTraversal.h
 * @date 08 Apr 2020
 * @author gratl
 *
 * adapted by henkel to use the sli_blk traversal
 */

#pragma once

#include <algorithm>

#include "LinkedCellTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/SlicedBlkBasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/C08CellHandler.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VerletListsCellsTraversal.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the multiple dimension sliced traversal.
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 */
template<class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class SlicedBlkTraversal : public SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                           public LinkedCellTraversalInterface<ParticleCell> {

public:
    /**
     * Constructor of the multiple dimension sliced traversal
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
     */
    explicit SlicedBlkTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                const double interactionLength, const std::array<double, 3> &cellLength) :
            SlicedBlkBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                        interactionLength,
                                                                                        cellLength),
            _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength,
                         this->_overlap) {}

    void traverseParticlePairs() override;

    DataLayoutOption getDataLayout() const override { return dataLayout; }

    bool getUseNewton3() const override { return useNewton3; }

    TraversalOption getTraversalType() const override { return TraversalOption::sli_blk; }

private:
    C08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;
};

template<class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void SlicedBlkTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
    auto &cells = *(this->_cells);
    this->slicedBlkTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
        auto id = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
        _cellHandler.processBaseCell(cells, id);
    });
}

}   // namespace autopas