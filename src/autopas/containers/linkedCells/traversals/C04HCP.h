/**
 * @file C04HCP.h
 * @author sabrinakrallmann
 * @date 30.03.2020
 */

#pragma once

#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/C08CellHandler.h"
#include "autopas/containers/linkedCells/traversals/LinkedCellTraversalInterface.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the c04 hcp traversal. The traversal is based on Tchipev, N. Algorithmic and Implementational
 * Optimizations of Molecular Dynamics Simulations for Process Engineering, Chapter 8 Outlook.
 *
 * The traversal uses the c04 base step performed on every single cell. Since
 * these steps overlap a domain coloring with four colors is applied. It differs from c04 in the shape of the colored
 * blocks. The chosen block-shape for the c08-base-steps in the c04HCP traversal is a 2x1x3-shape cuboid.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class C04HCP : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
               public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of c04hcp
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length.
   * @param cellLength cell length.
   */
  C04HCP(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor, const double interactionLength,
         const std::array<double, 3> &cellLength)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap),
        _end(utils::ArrayMath::subScalar(utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension), 1l)) {}

  void traverseParticlePairs() override;

  TraversalOption getTraversalType() const override { return TraversalOption::c04HCP; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

  bool isApplicable() const override {  // TODO: Individual Restrictions
    if (dataLayout == DataLayoutOption::cuda) {
      return false;
    }
    const double minLength = *std::min_element(this->_cellLength.cbegin(), this->_cellLength.cend());
    const unsigned long minDim = *std::min_element(this->_cellsPerDimension.cbegin(), this->_cellsPerDimension.cend());

    return minLength >= this->_interactionLength and minDim > 3;
  }

 private:
  void traverseSingleColor(std::vector<ParticleCell> &cells, int color);

  void processBasePack6(std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex);

  C08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;

  const std::array<long, 3> _end;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04HCP<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBasePack6(
    std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex) {
  using utils::ThreeDimensionalMapping::threeToOneD;
  std::array<long, 3> index{};
  const std::array<long, 3> signedDims = utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension);

  for (long z = 0; z < 3; ++z) {  // go through the six cells
    for (long x = 0; x < 2; ++x) {
      index[0] = base3DIndex[0] + x;
      index[1] = base3DIndex[1];
      index[2] = base3DIndex[2] + z;

      bool isIn = true;
      for (int d = 0; d < 3; ++d) {
        isIn &= (index[d] >= 0l) and (index[d] <= (_end[d] - this->_overlap[d]));  // prevent using overlapping cells
      }

      if (isIn) {  // skip cells outside radius
        const unsigned long ulIndex = threeToOneD(index, signedDims);
        _cellHandler.processBaseCell(cells, ulIndex);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04HCP<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
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

/**
 * Go through one color and search for cuboids belonging to the specified color.
 * Uses shifts to go through the different dimensions and prevent overlapping of the cuboids.
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 * @param cells
 * @param color
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04HCP<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseSingleColor(
    std::vector<ParticleCell> &cells, int color) {
  // determine a starting point of one of the grids
  std::array<long, 3> startOfThisColor{};  // coordinates: {x,y,z}

  // different starting points for different colors
  // some colors are starting outside the grid because only part of their cuboids are part of the grid
  // this way the starting points of sticking out cuboids can be determined as well
  switch (color) {
    case 0:
      startOfThisColor = {0l, 0l, 0l};
      break;
    case 1:
      startOfThisColor = {-4l, 0l, 1l};
      break;
    case 2:
      startOfThisColor = {-4l, 0l, -2l};
      break;
    case 3:
      startOfThisColor = {-2l, 0l, -1l};
      break;
  }

  // to fix compiler complaints about perfectly nested loop.
  long startX = startOfThisColor[0], endX = _end[0];
  long startY = startOfThisColor[1], endY = _end[1];
  long startZ = startOfThisColor[2], endZ = _end[2];

  // iterate over cartesian grid
#if defined(AUTOPAS_OPENMP)
#pragma omp for schedule(dynamic, 1) collapse(3) nowait
#endif
  for (long z = startZ; z < endZ; z += 4) {
    for (long y = startY; y < endY; y++) {
      for (long x = startX; x < (endX + 4);
           x += 6) {  // color starts every 6th column again, the +4 is needed to prevent ending too early, since it
                      // will be shifted back inside the loop
        long x_index = x;
        // shift on x-axis according to z-value
        // first: no shift
        // second: -4 shift
        // third: -2 shift
        // fourth: go back to first
        switch ((z - startZ) % 12 / 4) {
          case 0:
            break;
          case 1:
            x_index -= 4;
            break;
          case 2:
            x_index -= 2;
            break;
        }
        if ((y - startY) % 2 != 0) {  // shift x-axis every second y-row
          x_index += 3;
        }
        processBasePack6(cells, {x_index, y, z});
      }
    }
  }
}

}  // namespace autopas
