/**
 * @file C04Traversal.h
 * @author C.Menges, based on tchipevn (original source:
 * ls1-mardyn/src/particleContainer/LinkedCellTraversals/C04CellPairTraversal.h)
 * @date 15.06.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/C08BasedTraversal.h"
#include "autopas/containers/linkedCells/traversals/C08CellHandler.h"
#include "autopas/containers/linkedCells/traversals/LinkedCellTraversalInterface.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
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
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class C04Traversal : public C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                     public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c04 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength Interaction length.
   * @param cellLength cell length.
   */
  C04Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
               const double interactionLength, const std::array<double, 3> &cellLength)
      : C08BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                 interactionLength, cellLength),
        _cellHandler(pairwiseFunctor, this->_cellsPerDimension, interactionLength, cellLength, this->_overlap),
        _end(utils::ArrayMath::subScalar(utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension), 1l)) {
    computeOffsets32Pack();
  }

  void traverseParticlePairs() override;

  TraversalOption getTraversalType() const override { return TraversalOption::c04; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }

  bool getUseNewton3() const override { return useNewton3; }

  /**
   * C04 traversals are usable, if cellSizeFactor >= 1.0 and there are at least 3 cells for each dimension.
   * @return
   */
  bool isApplicable() const override {
    if (dataLayout == DataLayoutOption::cuda) {
      return false;
    }
    const double minLength = *std::min_element(this->_cellLength.cbegin(), this->_cellLength.cend());
    const unsigned long minDim = *std::min_element(this->_cellsPerDimension.cbegin(), this->_cellsPerDimension.cend());

    return minLength >= this->_interactionLength and minDim > 3;
  }

 private:
  void traverseSingleColor(std::vector<ParticleCell> &cells, int color);

  void processBasePack32(std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex);

  constexpr void computeOffsets32Pack();

  constexpr long parity(long x, long y, long z) const { return (x + y + z + 24) % 8; }

  std::array<std::array<long, 3>, 32> _cellOffsets32Pack;

  C08CellHandler<ParticleCell, PairwiseFunctor, dataLayout, useNewton3> _cellHandler;

  const std::array<long, 3> _end;
};

/**
 * Computes the barriers of the aggregation of cells for each color
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
constexpr void C04Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::computeOffsets32Pack() {
  using std::make_pair;
  using utils::ThreeDimensionalMapping::threeToOneD;

  unsigned int i = 0;
  long z = 0l;
  _cellOffsets32Pack[i++] = {1l, 1l, z};
  _cellOffsets32Pack[i++] = {1l, 2l, z};
  _cellOffsets32Pack[i++] = {2l, 1l, z};
  _cellOffsets32Pack[i++] = {2l, 2l, z};

  // z = 1ul; z = 2ul
  for (z = 1l; z < 3l; ++z) {
    for (long y = 0l; y < 4l; y++) {
      for (long x = 0l; x < 4l; x++) {
        if ((x == 0l and y == 0l) or (x == 3l and y == 0l) or (x == 0l and y == 3l) or (x == 3l and y == 3l)) {
          continue;
        }
        _cellOffsets32Pack[i++] = {x, y, z};
      }
    }
  }

  z = 3ul;
  _cellOffsets32Pack[i++] = {1l, 1l, z};
  _cellOffsets32Pack[i++] = {1l, 2l, z};
  _cellOffsets32Pack[i++] = {2l, 1l, z};
  _cellOffsets32Pack[i++] = {2l, 2l, z};

  if (i != 32) {
    utils::ExceptionHandler::exception("Internal error: Wrong number of offsets (expected: 32, actual: {})", i);
  }
}

/**
 * Goes through the cells aggregated by one color and processes the particles in each cell that is part of the
 * aggregation
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 * @param cells
 * @param base3DIndex
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processBasePack32(
    std::vector<ParticleCell> &cells, const std::array<long, 3> &base3DIndex) {
  using utils::ThreeDimensionalMapping::threeToOneD;
  std::array<long, 3> index;
  const std::array<long, 3> signedDims = utils::ArrayUtils::static_cast_array<long>(this->_cellsPerDimension);

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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void C04Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseSingleColor(
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

  // to fix compiler complaints about perfectly nested loop.
  const long startX = startOfThisColor[0], endX = _end[0];
  const long startY = startOfThisColor[1], endY = _end[1];
  const long startZ = startOfThisColor[2], endZ = _end[2];

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
      }
    }
  }
}

}  // namespace autopas
