/**
 * @file C01Traversal.h
 * @author nguyen
 * @date 16.09.2018
 */

#pragma once

#include "LinkedCellTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/C01BasedTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class provides the c01 traversal.
 *
 * The traversal uses the c01 base step performed on every single cell.
 * newton3 cannot be applied!
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class C01Traversal : public C01BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>,
                     public LinkedCellTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction (incl. halo).
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param cutoff cutoff radius
   * @param cellLength cell length in CellBlock3D
   */
  explicit C01Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                        const double cutoff = 1.0, const std::array<double, 3> &cellLength = {1.0, 1.0, 1.0})
      : C01BasedTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>(dims, pairwiseFunctor, cutoff,
                                                                                 cellLength),
        _cellFunctor(pairwiseFunctor, cutoff) {
    computeOffsets();
  }

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  /**
   * @copydoc LinkedCellTraversalInterface::traverseCellPairs()
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override;

  /**
   * C01 traversals are only usable if useNewton3 is disabled.
   *
   * This is because the cell functor in the c01 traversal is hardcoded to not allow newton 3 even if only one thread is
   * used.
   *
   * @return
   */
  bool isApplicable() override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    if (DataLayout == DataLayoutOption::cuda) {
      return (not useNewton3) and (nDevices > 0);
    } else {
      return not useNewton3;
    }
  }

  TraversalOption getTraversalType() override { return TraversalOption::c01; }

 private:
  /**
   * Computes all interactions between the base
   * cell and adjacent cells.
   * @param cells vector of all cells.
   * @param x X-index of base cell.
   * @param y Y-index of base cell.
   * @param z Z-index of base cell.
   */
  inline void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  /**
   * Pairs for processBaseCell().
   */
  std::vector<std::pair<long, std::array<double, 3>>> _cellOffsets;

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, false, false>
      _cellFunctor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::computeOffsets() {
  _cellOffsets.reserve(this->_overlap[0] * this->_overlap[1] * this->_overlap[2] * 3);

  const auto cutoffSquare(this->_cutoff * this->_cutoff);

  for (long z = -this->_overlap[0]; z <= static_cast<long>(this->_overlap[0]); ++z) {
    for (long y = -this->_overlap[1]; y <= static_cast<long>(this->_overlap[1]); ++y) {
      for (long x = -this->_overlap[2]; x <= static_cast<long>(this->_overlap[2]); ++x) {
        std::array<double, 3> pos = {};
        pos[0] = std::max(0l, (std::abs(x) - 1l)) * this->_cellLength[0];
        pos[1] = std::max(0l, (std::abs(y) - 1l)) * this->_cellLength[1];
        pos[2] = std::max(0l, (std::abs(z) - 1l)) * this->_cellLength[2];
        const double distSquare = ArrayMath::dot(pos, pos);
        if (distSquare <= cutoffSquare) {
          const long offset = (z * this->_cellsPerDimension[1] + y) * this->_cellsPerDimension[0] + x;
          _cellOffsets.push_back(std::make_pair(offset, ArrayMath::normalize(pos)));
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  unsigned long baseIndex = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];

  const size_t num_pairs = this->_cellOffsets.size();
  for (size_t j = 0; j < num_pairs; ++j) {
    const unsigned long otherIndex = baseIndex + this->_cellOffsets[j].first;
    ParticleCell &otherCell = cells[otherIndex];

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseCell);
    } else {
      std::array<double, 3> r = this->_cellOffsets[j].second;
      this->_cellFunctor.processCellPair(baseCell, otherCell, r);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
inline void C01Traversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells) {
  if (not this->isApplicable()) {
    utils::ExceptionHandler::exception(
        "The C01 traversal cannot work with enabled newton3 (unless only one thread is used)!");
  }
  this->c01Traversal([&](unsigned long x, unsigned long y, unsigned long z) { this->processBaseCell(cells, x, y, z); });
}

}  // namespace autopas
