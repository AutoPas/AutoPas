/**
 * @file C08BasedTraversal.h
 * @author F. Gratl
 * @date 5/4/18
 */

#pragma once

#include <utils/ThreeDimensionalMapping.h>
#include "CellPairTraversal.h"

namespace autopas {

template<class ParticleCell, class CellFunctor>
class C08BasedTraversal : public CellPairTraversals<ParticleCell, CellFunctor> {
 public:
  /**
   * Constructor of the c08 traversal.
   * @param cells The cells through which the traversal should traverse.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param cellfunctor The cell functor that defines the interaction of
   * particles between two different cells.
   */
  explicit C08BasedTraversal(std::vector<ParticleCell> &cells,
                        const std::array<unsigned long, 3> &dims,
                        CellFunctor *cellfunctor)
      : CellPairTraversals<ParticleCell, CellFunctor>(cells, dims,
                                                      cellfunctor) {
    computeOffsets();
  }

 protected:
  void processBaseCell(unsigned long baseIndex) const;
  void computeOffsets();

  std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets;
  std::array<unsigned long, 8> _cellOffsets;
};

template <class ParticleCell, class CellFunctor>
inline void C08BasedTraversal<ParticleCell, CellFunctor>::processBaseCell(
    unsigned long baseIndex) const {
  using std::pair;

  const int num_pairs = _cellPairOffsets.size();
  for (int j = 0; j < num_pairs; ++j) {
    pair<unsigned long, unsigned long> current_pair = _cellPairOffsets[j];

    unsigned long offset1 = current_pair.first;
    unsigned long cellIndex1 = baseIndex + offset1;

    unsigned long offset2 = current_pair.second;
    unsigned long cellIndex2 = baseIndex + offset2;

    ParticleCell &cell1 = this->_cells->at(cellIndex1);
    ParticleCell &cell2 = this->_cells->at(cellIndex2);

    if (cellIndex1 == cellIndex2) {
      this->_cellFunctor->processCell(cell1);
    } else {
      this->_cellFunctor->processCellPair(cell1, cell2);
    }
  }
}

template <class ParticleCell, class CellFunctor>
inline void C08BasedTraversal<ParticleCell, CellFunctor>::computeOffsets() {
  using ThreeDimensionalMapping::threeToOneD;
  using std::make_pair;

  unsigned long o =
      threeToOneD(0ul, 0ul, 0ul, this->_cellsPerDimension);  // origin
  unsigned long x = threeToOneD(
      1ul, 0ul, 0ul, this->_cellsPerDimension);  // displacement to the right
  unsigned long y =
      threeToOneD(0ul, 1ul, 0ul, this->_cellsPerDimension);  // displacement ...
  unsigned long z = threeToOneD(0ul, 0ul, 1ul, this->_cellsPerDimension);
  unsigned long xy = threeToOneD(1ul, 1ul, 0ul, this->_cellsPerDimension);
  unsigned long yz = threeToOneD(0ul, 1ul, 1ul, this->_cellsPerDimension);
  unsigned long xz = threeToOneD(1ul, 0ul, 1ul, this->_cellsPerDimension);
  unsigned long xyz = threeToOneD(1ul, 1ul, 1ul, this->_cellsPerDimension);

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