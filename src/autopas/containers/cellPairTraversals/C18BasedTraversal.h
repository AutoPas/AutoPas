/**
 * @file C18BasedTraversal.h
 * @author nguyen
 * @date 06.09.18
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/VerletListsTraversal.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {

/**
 * This class provides the base for traversals using the c18 base step.
 *
 * The base step processBaseCell() computes all interactions
 * between the base cell and adjacent cells with greater a ID.
 * After executing the base step on all cells all pairwise interactions for
 * all cells are done.
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam useSoA
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
class C18BasedTraversal : public CellPairTraversal<ParticleCell>,
                          public VerletListsTraversal<PairwiseFunctor, useNewton3> {
 public:
  /**
   * Constructor of the c18 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit C18BasedTraversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims),
        VerletListsTraversal<PairwiseFunctor, useNewton3>(pairwiseFunctor),
        _cellFunctor(
            CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3>(
                pairwiseFunctor)) {
    computeOffsets();
  }

 protected:
  /**
   * Computes all interactions between the base
   * cell and adjacent cells with greater a ID.
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
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, useSoA, useNewton3> _cellFunctor;

  /**
   * Pairs for processBaseCell(). 3x3 Array for each
   * special case in x and y direction.
   */
  std::array<std::array<std::vector<unsigned long>, 3>, 3> _cellOffsets;
};

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C18BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::processBaseCell(
    std::vector<ParticleCell> &cells, unsigned long x, unsigned long y, unsigned long z) {
  unsigned long baseIndex = ThreeDimensionalMapping::threeToOneD(x, y, z, this->_cellsPerDimension);

  unsigned int xArray;
  if (x == 0) {
    xArray = 0;
  } else if (x < this->_cellsPerDimension[0] - 1) {
    xArray = 1;
  } else {
    xArray = 2;
  }

  unsigned int yArray;
  if (y == 0) {
    yArray = 0;
  } else if (y < this->_cellsPerDimension[1] - 1) {
    yArray = 1;
  } else {
    yArray = 2;
  }

  ParticleCell &baseCell = cells[baseIndex];
  std::vector<unsigned long> &offsets = _cellOffsets[yArray][xArray];
  const int num_pairs = offsets.size();
  for (int j = 0; j < num_pairs; ++j) {
    unsigned long otherIndex = baseIndex + offsets[j];
    ParticleCell &otherCell = cells[otherIndex];

    if (baseIndex == otherIndex) {
      this->_cellFunctor.processCell(baseCell);
    } else {
      this->_cellFunctor.processCellPair(baseCell, otherCell);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, bool useSoA, bool useNewton3>
inline void C18BasedTraversal<ParticleCell, PairwiseFunctor, useSoA, useNewton3>::computeOffsets() {
  for (int z = 0; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        int offset = (z * this->_cellsPerDimension[1] + y) * this->_cellsPerDimension[0] + x;

        if (offset >= 0) {
          unsigned long uoffset = (unsigned long)offset;
          // add to each applicable special case
          for (int yArray = -1; yArray <= 1; ++yArray) {
            if (std::abs(yArray + y) <= 1) {
              for (int xArray = -1; xArray <= 1; ++xArray) {
                if (std::abs(xArray + x) <= 1) {
                  _cellOffsets[yArray + 1][xArray + 1].push_back(uoffset);
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace autopas
