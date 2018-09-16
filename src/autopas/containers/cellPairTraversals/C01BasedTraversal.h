/**
 * @file C01BasedTraversal.h
 * @author nguyen
 * @date 16.09.18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/VerletListsTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c01 base step.
 *
 * The base step processBaseCell() computes all interactions
 * between the base cell and adjacent cells.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA>
class C01BasedTraversal : public CellPairTraversal<ParticleCell>, public VerletListsTraversal<PairwiseFunctor, false> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C01BasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims),
	VerletListsTraversal<PairwiseFunctor, false> (pairwiseFunctor),
        _cellFunctor(CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, false, false>(pairwiseFunctor)) {
    computeOffsets();
  }

 protected:
  /**
   * Computes all interactions between the base
   * cell and adjacent cells.
   * @param cells vector of all cells.
   * @param x x of base cell
   * @param y y of base cell
   * @param z z of base cell
   */
  void processBaseCell(std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z);

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, false, false> _cellFunctor;

  /**
   * Pairs for processBaseCell().
   */
  std::vector<int> _cellOffsets;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  unsigned long baseIndex = ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);
  ParticleCell &baseCell = cells[baseIndex];

  const int num_pairs = _cellOffsets.size();
  for (int j = 0; j < num_pairs; ++j) {
    unsigned long otherIndex = baseIndex + _cellOffsets[j];
    ParticleCell &otherCell = cells[otherIndex];

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseCell);
    } else {
      this->_cellFunctor.processCellPair(baseCell, otherCell);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA>
inline void C01BasedTraversal<ParticleCell, PairwiseFunctor, useSoA>::computeOffsets() {
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        int offset = (z * this->_cellsPerDimension[1] + y) * this->_cellsPerDimension[0] + x;
	_cellOffsets.push_back(offset);
      }
    }
  }
}

}  // namespace autopas
