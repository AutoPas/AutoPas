/**
 * @file C18Traversal.h
 * @author nguyen
 * @date 06.09.2018
 */

#pragma once

#include "LinkedCellTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C18BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c18 traversal.
 *
 * The traversal uses the c18 base step performed on every single cell. Since
 * these steps overlap a domain coloring with eighteen colors is applied.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C18Traversal : public C18BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>,
                     public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff Cutoff radius.
   * @param cellLength cell length.
   */
  explicit C18Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                        const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : C18BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor, cutoff,
                                                                                 cellLength),
        _cellFunctor(pairwiseFunctor, cutoff) {
    computeOffsets();
  }

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * Computes all interactions between the base
   * cell and adjacent cells with greater a ID.
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  TraversalOption getTraversalType() override { return TraversalOption::c18; }

  /**
   * C18 traversal is always usable.
   * @return
   */
  bool isApplicable() override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    if (DataLayout == DataLayoutOption::cuda)
      return nDevices > 0;
    else
      return true;
  }

 private:
  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3> _cellFunctor;

  /**
   * Pairs for processBaseCell(). overlap[0] x overlap[1] Array for each special case in x and y direction.
   */
  std::vector<std::vector<std::vector<std::pair<unsigned long, std::array<double, 3>>>>> _cellOffsets;

  /**
   * Returns the index in the offsets array for the given position.
   * @param pos current position in dimension dim
   * @param dim current dimension
   * @return Index for the _cellOffsets Array.
   */
  unsigned int getIndex(unsigned long pos, unsigned int dim) const;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C18Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::computeOffsets() {
  _cellOffsets.resize(
      2 * this->_overlap[1] + 1,
      std::vector<std::vector<std::pair<unsigned long, std::array<double, 3>>>>(2 * this->_overlap[0] + 1));
  const std::array<long, 3> _overlap_s = {static_cast<long>(this->_overlap[0]), static_cast<long>(this->_overlap[1]),
                                          static_cast<long>(this->_overlap[2])};

  const auto cutoffSquare(this->_cutoff * this->_cutoff);

  for (long z = 0l; z <= _overlap_s[2]; ++z) {
    for (long y = -_overlap_s[1]; y <= _overlap_s[1]; ++y) {
      for (long x = -_overlap_s[0]; x <= _overlap_s[0]; ++x) {
        const long offset = (z * this->_cellsPerDimension[1] + y) * this->_cellsPerDimension[0] + x;
        if (offset >= 0l) {
          // add to each applicable special case
          for (long yArray = -_overlap_s[1]; yArray <= _overlap_s[1]; ++yArray) {
            if (std::abs(yArray + y) <= _overlap_s[1]) {
              for (long xArray = -_overlap_s[0]; xArray <= _overlap_s[0]; ++xArray) {
                if (std::abs(xArray + x) <= _overlap_s[0]) {
                  std::array<double, 3> pos = {};
                  pos[0] = std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0];
                  pos[1] = std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1];
                  pos[2] = std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2];
                  // calculate distance between base cell and other cell
                  const double distSquare = ArrayMath::dot(pos, pos);
                  // only add cell offset if cell is within cutoff radius
                  if (distSquare <= cutoffSquare) {
                    auto uoffset = static_cast<unsigned long>(offset);
                    _cellOffsets[yArray + _overlap_s[1]][xArray + _overlap_s[0]].push_back(
                        std::make_pair(offset, ArrayMath::normalize(pos)));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
unsigned int C18Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::getIndex(unsigned long pos,
                                                                                           unsigned int dim) const {
  unsigned long index;
  if (pos < this->_overlap[dim]) {
    index = pos;
  } else if (pos < this->_cellsPerDimension[dim] - this->_overlap[dim]) {
    index = this->_overlap[dim];
  } else {
    index = pos - this->_cellsPerDimension[dim] + 2 * this->_overlap[dim] + 1ul;
  }
  return index;
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void C18Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  const unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  const unsigned int xArray = getIndex(x, 0);

  const unsigned int yArray = getIndex(y, 1);

  ParticleCell &baseCell = cells[baseIndex];
  std::vector<std::pair<unsigned long, std::array<double, 3>>> &offsets = this->_cellOffsets[yArray][xArray];
  const size_t num_pairs = offsets.size();
  for (size_t j = 0; j < num_pairs; ++j) {
    unsigned long otherIndex = baseIndex + offsets[j].first;
    ParticleCell &otherCell = cells[otherIndex];

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseCell);
    } else {
      this->_cellFunctor.processCellPair(baseCell, otherCell, offsets[j].second);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C18Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  this->c18Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
