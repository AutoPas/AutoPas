/*
 * SlicedTraversal.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#pragma once

#include <utils/WrapOpenMP.h>
#include <algorithm>
#include "CellPairTraversal.h"
#include "utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the sliced traversal.
 * @todo enhance documentation
 * @tparam ParticleCell the type of cells
 * @tparam CellFunctor the cell functor that defines the interaction of the
 * particles of two specific cells
 */
template <class ParticleCell, class CellFunctor>
class SlicedTraversal : public CellPairTraversals<ParticleCell, CellFunctor> {
 public:
  /**
   * constructor of the sliced traversal
   * @param cells the cells through which the traversal should traverse
   * @param dims the dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction
   * @param cellfunctor the cell functor that defines the interaction of
   * particles between two different cells
   */
  explicit SlicedTraversal(std::vector<ParticleCell> &cells,
                           const std::array<unsigned long, 3> &dims,
                           CellFunctor *cellfunctor)
      : CellPairTraversals<ParticleCell, CellFunctor>(cells, dims,
                                                      cellfunctor) {
    computeOffsets();
  }
  // documentation in base class
  void traverseCellPairs() override;

 private:
  void processBaseCell(unsigned long baseIndex) const;
  void computeOffsets();

  std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets;
  std::array<unsigned long, 8> _cellOffsets;
};

template <class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::processBaseCell(
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
inline void SlicedTraversal<ParticleCell, CellFunctor>::computeOffsets() {
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

template <class ParticleCell, class CellFunctor>
inline void SlicedTraversal<ParticleCell, CellFunctor>::traverseCellPairs() {
  //  using std::array;
  //  array<unsigned long, 3> endid;
  //  for (int d = 0; d < 3; ++d) {
  //    endid[d] = this->_cellsPerDimension[d] - 1;
  //  }
  //
  //  for (unsigned long z = 0; z < endid[2]; ++z) {
  //    for (unsigned long y = 0; y < endid[1]; ++y) {
  //      for (unsigned long x = 0; x < endid[0]; ++x) {
  //        unsigned long ind =
  //            ThreeDimensionalMapping::threeToOneD(x, y, z,
  //            this->_cellsPerDimension);
  //        processBaseCell(ind);
  //      }
  //    }
  //  }

  // 0) TODO: check if applicable

  // 1) split domain across its longest dimension

  // find dimension with most cells

  std::vector<std::pair<int, unsigned long>> mapDimLength{
      {0, this->_cellsPerDimension[0]},
      {1, this->_cellsPerDimension[1]},
      {2, this->_cellsPerDimension[2]},
  };

  // sorts longest dimension to index 0
  std::sort(
      mapDimLength.begin(), mapDimLength.end(),
      [](std::pair<int, unsigned long> a, std::pair<int, unsigned long> b) {
        return a.second > b.second;
      });

  auto numSlices = autopas_get_num_threads();

  unsigned long sliceThickness = mapDimLength[0].second / numSlices;

  unsigned long cellsPerSlice = mapDimLength[1].second * mapDimLength[2].second;

  std::vector<autopas_lock_t *> locks;
  locks.reserve(numSlices);

  for (auto i = 0; i < numSlices - 1; ++i) {
    autopas_init_lock(locks[i]);
  }

#if defined(_OPENMP)
#pragma omp parallel for schedule(static, 1)
#endif
  for (auto slice = 0; slice < numSlices; ++slice) {
    // lock the starting layer
    if (slice > 0) autopas_set_lock(locks[slice - 1]);
    for (auto dimSlice = 0; dimSlice < sliceThickness; ++dimSlice) {
      // at the last layer request next lock the last layer has no next
      if (slice != numSlices - 1 && dimSlice == sliceThickness - 1) autopas_set_lock(locks[slice]);
      for (auto dimMedium = 0; dimMedium < mapDimLength[1].second;
           ++dimMedium) {
        for (auto dimShort = 0; dimShort < mapDimLength[2].second; ++dimShort) {
          processBaseCell(slice * sliceThickness + dimMedium * mapDimLength[2].second + dimShort);
        }
      }
      // at the end of the first layer release the lock
      if (slice > 0 && dimSlice == 0)
        autopas_unset_lock(locks[slice - 1]);
      else if (slice != numSlices - 1 &&dimSlice == sliceThickness - 1)
        autopas_unset_lock(locks[slice]);
    }
  }
}

} /* namespace autopas */
