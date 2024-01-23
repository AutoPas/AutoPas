/**
* @file LCC04Traversal3B.h
* @date 21 Nov 2023
* @author N. Deng
*/

#pragma once

#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/containers/cellTraversals/C08BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/LCC08CellHandler3B.h"
#include "autopas/containers/linkedCells/traversals/LCTraversalInterface.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
* This class provides the c04 traversal.
*
* The traversal uses the c04 base step performed on every single cell. Since
* these steps overlap a domain coloring with four colors is applied.
*
* @tparam ParticleCell the type of cells
* @tparam Functor The functor that defines the interaction of three particles.
* @tparam useSoA
* @tparam useNewton3
*/
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC04Traversal3B
   : public C08BasedTraversal<ParticleCell, Functor, InteractionTypeOption::threeBody, dataLayout, useNewton3>,
     public LCTraversalInterface<ParticleCell> {
public:
 /**
  * Constructor of the c04 traversal.
  * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
  * y and z direction.
  * @param Functor The functor that defines the interaction of two particles.
  * @param interactionLength Interaction length.
  * @param cellLength cell length.
  */
 LCC04Traversal3B(const std::array<unsigned long, 3> &dims, Functor *functor,
                const double interactionLength, const std::array<double, 3> &cellLength)
     : C08BasedTraversal<ParticleCell, Functor, InteractionTypeOption::threeBody, dataLayout, useNewton3>(
           dims, functor, interactionLength, cellLength),
       _cellOffsets32Pack(computeOffsets32Pack()),
       _cellHandler(functor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap),
       _end(utils::ArrayMath::subScalar(utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension),
                                        1l)) {}

 void traverseParticleTriplets() override;

 [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c04_3b; }

 [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

 [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

 /**
  * C04 traversals are usable, if cellSizeFactor >= 1.0 and there are at least 3 cells for each dimension.
  * @return information about applicability
  */
 [[nodiscard]] bool isApplicable() const override {
   // The cellsize cannot be smaller than the cutoff, if OpenMP is used.
   // Also see: https://github.com/AutoPas/AutoPas/issues/464
   const double minLength = *std::min_element(this->_cellLength.cbegin(), this->_cellLength.cend());
   const unsigned long minDim = *std::min_element(this->_cellsPerDimension.cbegin(), this->_cellsPerDimension.cend());

   return minLength >= this->_interactionLength and minDim > 3;
 }

 /**
  * @copydoc autopas::CellTraversal::setSortingThreshold()
  */
 void setSortingThreshold(size_t sortingThreshold) override { _cellHandler.setSortingThreshold(sortingThreshold); }

private:
 void traverseSingleColor(std::vector<ParticleCell> &cells, int color);

 void processBasePack32(std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex);

 constexpr auto computeOffsets32Pack() const;

 [[nodiscard]] constexpr long parity(long x, long y, long z) const { return (x + y + z + 24) % 8; }

 std::array<std::array<long, 3>, 32> _cellOffsets32Pack;

 LCC08CellHandler3B<ParticleCell, Functor, dataLayout, useNewton3> _cellHandler;

 const std::array<long, 3> _end;
};

/**
* Computes the barriers of the aggregation of cells for each color
*
* @tparam ParticleCell
* @tparam Functor
* @tparam dataLayout
* @tparam useNewton3
*/
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
constexpr auto LCC04Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::computeOffsets32Pack() const {
 using std::make_pair;
 using utils::ThreeDimensionalMapping::threeToOneD;

 std::array<std::array<long, 3>, 32> cellOffsets32Pack = {};

 unsigned int i = 0;
 long z = 0l;
 cellOffsets32Pack[i++] = {1l, 1l, z};
 cellOffsets32Pack[i++] = {1l, 2l, z};
 cellOffsets32Pack[i++] = {2l, 1l, z};
 cellOffsets32Pack[i++] = {2l, 2l, z};

 // z = 1ul; z = 2ul
 for (z = 1l; z < 3l; ++z) {
   for (long y = 0l; y < 4l; y++) {
     for (long x = 0l; x < 4l; x++) {
       if ((x == 0l and y == 0l) or (x == 3l and y == 0l) or (x == 0l and y == 3l) or (x == 3l and y == 3l)) {
         continue;
       }
       cellOffsets32Pack[i++] = {x, y, z};
     }
   }
 }

 z = 3ul;
 cellOffsets32Pack[i++] = {1l, 1l, z};
 cellOffsets32Pack[i++] = {1l, 2l, z};
 cellOffsets32Pack[i++] = {2l, 1l, z};
 cellOffsets32Pack[i++] = {2l, 2l, z};

 if (i != 32) {
   utils::ExceptionHandler::exception("Internal error: Wrong number of offsets (expected: 32, actual: {})", i);
 }

 return cellOffsets32Pack;
}

/**
* Goes through the cells aggregated by one color and processes the particles in each cell that is part of the
* aggregation by using the barriers saved in _cellOffset32Pack.
*
* @tparam ParticleCell
* @tparam Functor
* @tparam dataLayout
* @tparam useNewton3
* @param cells
* @param base3DIndex
*/
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
void LCC04Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::processBasePack32(
   std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex) {
 using utils::ThreeDimensionalMapping::threeToOneD;
 std::array<long, 3> index;
 const std::array<long, 3> signedDims = utils::ArrayUtils::static_cast_copy_array<long>(this->_cellsPerDimension);

 for (auto Offset32Pack : _cellOffsets32Pack) {
   // compute 3D index
   bool isIn = true;
   for (int d = 0; d < 3; ++d) {
     index[d] = base3DIndex[d] + Offset32Pack[d];
     isIn &= (index[d] >= 0l) and (index[d] < _end[d]);
   }

   if (isIn) {
     const unsigned long ulIndex = threeToOneD(index, signedDims);
     _cellHandler.processBaseCell(cells, ulIndex);
   }
 }
}

/**
*  Go through one color and search for blocks belonging to the specified color.
*  Uses two cartesian grids that are overlapping gridwise but not blockwise.
*
* @tparam ParticleCell
* @tparam Functor
* @tparam dataLayout
* @tparam useNewton3
*/
template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
void LCC04Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::traverseParticleTriplets() {
 auto &cells = *(this->_cells);
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
 {
   for (int color = 0; color < 4; ++color) {
     traverseSingleColor(cells, color);

#if defined(AUTOPAS_OPENMP)
     if (color < 3) {
#pragma omp barrier
     }
#endif
   }
 }  // close parallel region
}

template <class ParticleCell, class Functor, DataLayoutOption::Value dataLayout, bool useNewton3>
void LCC04Traversal3B<ParticleCell, Functor, dataLayout, useNewton3>::traverseSingleColor(
   std::vector<ParticleCell> &cells, int color) {
 // we need to traverse one body-centered cubic (BCC) grid, which consists of two cartesian grids

 // colors 0 and 2 form one cartesian grid
 // colors 1 and 3 form another cartesian grid, whose origin is shifted by (2,2,2)

 // determine a starting point of one of the grids
 std::array<long, 3> startOfThisColor{};

 switch (color % 2) {
   case 0:
     // colours 0 and 2
     startOfThisColor = {-2l, -2l, -2l};
     break;
   case 1:
     // colours 1 and 3
     startOfThisColor = {0l, 0l, 0l};
     break;
 }

 // calculate whether the calculated starting point is part of the color
 long correctParity = parity(startOfThisColor[0], startOfThisColor[1], startOfThisColor[2]);
 if (color >= 2) {
   correctParity += 4;
 }

 // to fix intel64 icpc compiler complaints about perfectly nested loop (tested with version 19.0.4.20190416).
 const long startX = startOfThisColor[0], endX = _end[0];
 const long startY = startOfThisColor[1], endY = _end[1];
 const long startZ = startOfThisColor[2], endZ = _end[2];

 auto counter = 0;
// first cartesian grid
// grids are interlinked: one grid fills the gaps in the other grid
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3) nowait
#endif
 for (long z = startZ; z < endZ; z += 4) {
   for (long y = startY; y < endY; y += 4) {
     for (long x = startX; x < endX; x += 4) {
       const long par = parity(x, y, z);

       if (par != correctParity) {
         continue;
       }

       const std::array<long, 3> base3DIndex = {x, y, z};
       processBasePack32(cells, base3DIndex);
       counter++;
     }
   }
 }
 std::cout << "C04 32 packs counter : " << counter << std::endl;
}

}  // namespace autopas
