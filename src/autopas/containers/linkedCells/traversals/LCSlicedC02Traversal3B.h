/**
* @file LCSlicedC02Traversal3B.h
*
* @date 22 Nov 2023
* @author N. Deng
*/

#pragma once

#include <algorithm>

#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/SlicedC02BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler3B.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
* This class provides the colored sliced traversal.
*
* The traversal finds the longest dimension of the simulation domain and cuts
* the domain into as many slices as possible along this dimension. Unlike the regular
* sliced traversal, this version uses a 2-coloring to prevent race conditions, instead of
* locking the starting layers.
*
* @tparam ParticleCell the type of cells
* @tparam Functor The functor that defines the interaction of three particles.
* @tparam DataLayout
* @tparam useNewton3
*/
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCSlicedC02Traversal3B
    : public SlicedC02BasedTraversal<ParticleCell, Functor, InteractionTypeOption::threeBody, dataLayout, useNewton3, true>,
     public LCTraversalInterface<ParticleCell> {
public:
 /**
  * Constructor of the colored sliced traversal.
  * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
  * y and z direction.
  * @param functor The functor that defines the interaction of three particles.
  * @param interactionLength Interaction length (cutoff + skin).
  * @param cellLength cell length.
  */
 explicit LCSlicedC02Traversal3B(const std::array<unsigned long, 3> &dims, Functor *functor,
                               const double interactionLength, const std::array<double, 3> &cellLength)
     : SlicedC02BasedTraversal<ParticleCell, Functor, InteractionTypeOption::threeBody, dataLayout, useNewton3, true>(
           dims, functor, interactionLength, cellLength),
       _cellHandler(functor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

 void traverseParticleTriplets() override;

 [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

 [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

 [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_sliced_c02_3b; }

 /**
  * @copydoc autopas::CellTraversal::setSortingThreshold()
  */
 void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

private:
 LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCSlicedC02Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::traverseParticleTriplets() {
 auto &cells = *(this->_cells);
 this->cSlicedTraversal([&](unsigned long x, unsigned long y, unsigned long z) {
   auto id = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
   _cellHandler.processBaseCell(cells, id);
 });
}

}  // namespace autopas
