/**
* @file LCC08Traversal3BOptimized.h
* @author N. Deng
* @date 28.10.2023
*/

#pragma once

#include "LCC08CellHandler3BOptimized.h"
#include "LCTraversalInterface.h"
#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
* This class provides the lc_c08 traversal for 3-body interactions.
*
* The traversal uses the c08 base step performed on every single cell.
* \image html C08.png "C08 base step in 2D. (dark blue cell = base cell)"
* Since these steps overlap a domain coloring with eight colors is applied.
* \image html C08_domain.png "C08 domain coloring in 2D. 4 colors are required."
*
* @tparam ParticleCell the type of cells
* @tparam functor The functor that defines the interaction of three particles.
* @tparam dataLayout
* @tparam useNewton3
*/
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC08Traversal3BOptimized
    : public C08BasedTraversal<ParticleCell, Functor, InteractionTypeOption::threeBody, dataLayout, useNewton3>,
      public LCTraversalInterface<ParticleCell> {
public:
 /**
  * Constructor of the lc_c08_3b traversal.
  * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
  * y and z direction.
  * @param functor The functor that defines the interaction of three particles.
  * @param interactionLength Interaction length (cutoff + skin).
  * @param cellLength cell length.
  */
 explicit LCC08Traversal3BOptimized(const std::array<unsigned long, 3> &dims, Functor *functor, const double interactionLength,
                           const std::array<double, 3> &cellLength)
     : C08BasedTraversal<ParticleCell, Functor, InteractionTypeOption::threeBody, dataLayout, useNewton3>(
           dims, functor, interactionLength, cellLength),
       _cellHandler(functor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap) {}

 void traverseParticleTriplets() override;

 [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c08_3b; }

 [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

 [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

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
 LCC08CellHandler3BOptimized<ParticleCell, Functor, dataLayout, useNewton3> _cellHandler;
};

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC08Traversal3BOptimized<ParticleCell, Functor, dataLayout, useNewton3>::traverseParticleTriplets() {
 auto &cells = *(this->_cells);
 this->c08Traversal([&](unsigned long x, unsigned long y, unsigned long z) {
   unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
   _cellHandler.processBaseCell(cells, baseIndex);
 });
}

}  // namespace autopas
