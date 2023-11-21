/**
 * @file LCSlicedTraversal3B.h
 *
 * @date 21 Nov 2023
 * @author N. Deng
 */

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/SlicedLockBasedTraversal.h"
#include "autopa"
//#include "autopas/containers/verletListsCellBased/verletListsCells/traversals/VLCTraversalInterface.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the (locked) sliced traversal.
 *
 * The traversal finds the longest dimension of the simulation domain and cuts
 * the domain into multiple slices along this dimension. Slices are
 * assigned to the threads in a round robin fashion. Each thread locks the cells
 * on the boundary wall to the previous slice with one lock. This lock is lifted
 * as soon the boundary wall is fully processed.
 *
 * @tparam ParticleCell the type of cells
 * @tparam Functor The functor that defines the interaction of three particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
    template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
    class LCSlicedTraversal3B : public SlicedLockBasedTraversal<ParticleCell, Functor, dataLayout, useNewton3, true>,
                              public LCTraversalInterface<ParticleCell> {
    public:
        /**
         * Constructor of the sliced traversal.
         * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
         * y and z direction.
         * @param functor The functor that defines the interaction of three particles.
         * @param interactionLength Interaction length (cutoff + skin).
         * @param cellLength cell length.
         */
        explicit LCSlicedTraversal3B(const std::array<unsigned long, 3> &dims, Functor *functor,
                                   const double interactionLength, const std::array<double, 3> &cellLength)
                : SlicedLockBasedTraversal<ParticleCell, Functor, dataLayout, useNewton3, true>(
                dims, functor, interactionLength, cellLength),
                  _cellHandler(functor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

        void traverseParticleTriplets() override;

        [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

        [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

        [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_sliced_3b; }


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
        LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3> _cellHandler;
    };

    template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
    inline void LCSlicedTraversal3B<ParticleCell, Functor, dataLayout, useNewton3>::traverseParticleTriplets() {
        auto &cells = *(this->_cells);
        this->slicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
            unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
            _cellHandler.processBaseCell(cells, baseIndex);
        });
    }

}  // namespace autopas
