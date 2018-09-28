/**
 * @file C08BasedTraversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c08 base step.
 *
 * The base step processBaseCell() computes one set of pairwise interactions
 * between two cells for each spatial direction based on the baseIndex.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C08BasedTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C08BasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims),
        _cellFunctor(
            CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3>(
                pairwiseFunctor)) {
    computeOffsets();
  }

 protected:
  /**
   * Computes one interaction for each spacial direction based on the lower left
   * frontal corner (=base index) of a 2x2x2 block of cells.
   * @param cells vector of all cells.
   * @param baseIndex Index respective to which box is constructed.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long baseIndex);
  /**
   * Computes pairs for the block used in processBaseCell()
   */
  void computeOffsets();

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3> _cellFunctor;

  /**
   * Pair sets for processBaseCell().
   */
  std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets;
  /**
   * Offsets for single cells
   */
  std::array<unsigned long, 8> _cellOffsets;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long baseIndex) {
  using std::pair;

  const int num_pairs = _cellPairOffsets.size();
  for (int j = 0; j < num_pairs; ++j) {
    pair<unsigned long, unsigned long> current_pair = _cellPairOffsets[j];

    unsigned long offset1 = current_pair.first;
    unsigned long cellIndex1 = baseIndex + offset1;

    unsigned long offset2 = current_pair.second;
    unsigned long cellIndex2 = baseIndex + offset2;

    ParticleCell &cell1 = cells[cellIndex1];
    ParticleCell &cell2 = cells[cellIndex2];

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor.processCell(cell1);
    } else {
      this->_cellFunctor.processCellPair(cell1, cell2);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C08BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::computeOffsets() {
  using std::make_pair;

  // origin
  unsigned long o = utils::ThreeDimensionalMapping::threeToOneD(0ul, 0ul, 0ul, this->_cellsPerDimension);
  // displacement to right ...
  unsigned long x = utils::ThreeDimensionalMapping::threeToOneD(1ul, 0ul, 0ul, this->_cellsPerDimension);
  unsigned long y = utils::ThreeDimensionalMapping::threeToOneD(0ul, 1ul, 0ul, this->_cellsPerDimension);
  unsigned long z = utils::ThreeDimensionalMapping::threeToOneD(0ul, 0ul, 1ul, this->_cellsPerDimension);
  unsigned long xy = utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 0ul, this->_cellsPerDimension);
  unsigned long yz = utils::ThreeDimensionalMapping::threeToOneD(0ul, 1ul, 1ul, this->_cellsPerDimension);
  unsigned long xz = utils::ThreeDimensionalMapping::threeToOneD(1ul, 0ul, 1ul, this->_cellsPerDimension);
  unsigned long xyz = utils::ThreeDimensionalMapping::threeToOneD(1ul, 1ul, 1ul, this->_cellsPerDimension);

  int i = 0;
  // if incrementing along X, the following order will be more cache-efficient:
  _cellPairOffsets[i++] = make_pair(o, o);
  _cellPairOffsets[i++] = make_pair(o, y);
  _cellPairOffsets[i++] = make_pair(y, z);
  _cellPairOffsets[i++] = make_pair(o, z);
  _cellPairOffsets[i++] = make_pair(o, yz);

  _cellPairOffsets[i++] = make_pair(x, yz);
  _cellPairOffsets[i++] = make_pair(x, y);
  _cellPairOffsets[i++] = make_pair(x, z);
  _cellPairOffsets[i++] = make_pair(o, x);
  _cellPairOffsets[i++] = make_pair(o, xy);
  _cellPairOffsets[i++] = make_pair(xy, z);
  _cellPairOffsets[i++] = make_pair(y, xz);
  _cellPairOffsets[i++] = make_pair(o, xz);
  _cellPairOffsets[i++] = make_pair(o, xyz);

  i = 0;
  _cellOffsets[i++] = o;
  _cellOffsets[i++] = y;
  _cellOffsets[i++] = z;
  _cellOffsets[i++] = yz;

  _cellOffsets[i++] = x;
  _cellOffsets[i++] = xy;
  _cellOffsets[i++] = xz;
  _cellOffsets[i++] = xyz;
}

}  // namespace autopas